//
//  FiniteDifference.cpp
//  FD
//
//  Created by Nicolas Buchwalder on 07.12.22.
//

#include <iostream>
#include "MonteCarlo.hpp"
#include "Option.hpp"
#include "Decomposer.hpp"
#include "FiniteDifference.hpp"
#include "IterativeSolver.hpp"
#include "LinearSolver.hpp"
#include "RandomNumberGenerator.hpp"

MonteCarlo::MonteCarlo(Option option, RandomNumberGeneratorEnum rng)
{
    std::tie(_ex, _payoff, _type, _S, _K, _T, _sigma, _r, _q, _divs, _add_params) = option.get_params();
    _option = option;
    _rng = rng;
}

MonteCarlo::MonteCarlo(Option option, RandomNumberGeneratorEnum rng, Model model, std::vector<double> Heston)
{
    std::tie(_ex, _payoff, _type, _S, _K, _T, _sigma, _r, _q, _divs, _add_params) = option.get_params();
    _option = option;
    _rng = rng;
    _Model = model;
    _Heston = Heston;
}

RandomNumberGenerator *MonteCarlo::getRandomNumberGenerator()
{
    switch (_rng)
    {
    case RandomNumberGeneratorEnum::LinearCongruentialRNG:
        return new LinearCongruential();

    case RandomNumberGeneratorEnum::InverseTransformRNG:
        return new InverseTransform();

    case RandomNumberGeneratorEnum::AcceptanceRejectionRNG:
        return new AcceptanceRejection();

    case RandomNumberGeneratorEnum::BoxMullerRNG:
        return new BoxMuller();

    default:
        return new BoxMuller();
    }
}

void MonteCarlo::validateParameters()
{
    if (!(_Model == Model::SDE || _Model == Model::Heston))
    {
        throw std::invalid_argument("Unknown model provided");
    }

    if (_Model == Model::Heston && _Heston.size() != 4)
    {
        throw std::invalid_argument("Model provided is Heston but Heston parameters are not defined");
    }
}

/**
 * Generates one MC path.
 */
double *MonteCarlo::generatePathSDE(double *randoms, std::size_t timeSteps)
{
    double *Sj = new double[timeSteps + 1];
    double dt = _T / timeSteps;
    Sj[0] = _S;
    for (std::size_t i = 1; i <= timeSteps; i++)
    {
        Sj[i] = Sj[i - 1] * exp((_r - _q - pow(_sigma, 2) / 2) * dt + _sigma * sqrt(dt) * randoms[i]);
    }

    return Sj;
}

double *MonteCarlo::generatePathHeston(double *randomSpots, double *randomVols, std::size_t timeSteps)
{
    // Unpack Heston parameters;
    double lambda = _Heston[0];
    double sigmalongterm = _Heston[1];
    double eeta = _Heston[2];
    double rho = _Heston[3];
    double Vlongterm = pow(sigmalongterm, 2);

    double *Sj = new double[timeSteps + 1];
    double *Vj = new double[timeSteps + 1];
    double dt = _T / timeSteps;
    Sj[0] = _S;
    Vj[0] = pow(_sigma, 2);

    for (std::size_t i = 1; i <= timeSteps; i++)
    {
        Sj[i] = Sj[i - 1] * exp((_r - _q - std::max(Vj[i - 1], 0.0) / 2) * dt + sqrt(std::max(Vj[i - 1], 0.0)) * sqrt(dt) * randomSpots[i]);
        Vj[i] = std::max(Vj[i - 1], 0.0) - lambda * (std::max(Vj[i - 1], 0.0) - Vlongterm) * dt + eeta * sqrt(std::max(Vj[i - 1], 0.0)) * sqrt(dt) * (rho * randomSpots[i] + sqrt(1 - pow(rho, 2)) * randomVols[i]);
    }

    return Sj;
}

std::vector<double> MonteCarlo::price_option(std::size_t numPaths, std::size_t timeSteps, bool include_greeks)
{
    MonteCarlo::validateParameters();
    bool isBarrierOption = _option.isBarrierOption();
    bool shouldBarrierBeHit = _type == OptionType::downin || _type == OptionType::upin;

    std::vector<double> price_return;

    std::size_t len = numPaths * timeSteps + 1;
    double *randoms = new double[len];
    RandomNumberGenerator *rng = getRandomNumberGenerator();
    for (std::size_t i = 0; i < len; i++)
    {
        randoms[i] = (*rng)();
    }

    double *path;
    double Vcap = 0;
    for (std::size_t i = 0; i < numPaths; i++)
    {
        path = generatePathSDE(&randoms[i * timeSteps], timeSteps);
        if (!isBarrierOption || (shouldBarrierBeHit && isBarrierHit(path, timeSteps)) // downin and upin options - Barrier should be hit and it is hit
            || (!shouldBarrierBeHit && !isBarrierHit(path, timeSteps))                // downout and upout options - Barrier should not be hit and it is not hit
        )
        {
            Vcap += exp(-_r * _T) * getPathPayoff(path[timeSteps]);
        }
    }

    Vcap = Vcap / numPaths;
    price_return.push_back(Vcap);
    return price_return;
}

std::vector<double> MonteCarlo::price_option_controlvariate(std::size_t numPaths, std::size_t timeSteps, bool include_greeks)
{
    MonteCarlo::validateParameters();
    bool isBarrierOption = _option.isBarrierOption();
    bool shouldBarrierBeHit = _type == OptionType::downin || _type == OptionType::upin;

    std::vector<double> price_return;

    std::size_t len = numPaths * timeSteps + 1;
    double *randoms = new double[len];
    RandomNumberGenerator *rng = getRandomNumberGenerator();
    for (std::size_t i = 0; i < len; i++)
    {
        randoms[i] = (*rng)();
    }

    double *path;
    double *S_arr = new double[numPaths];
    double *V_arr = new double[numPaths];
    double Scap = 0;
    double Vcap = 0;

    for (std::size_t i = 0; i < numPaths; i++)
    {
        path = generatePathSDE(&randoms[i * timeSteps], timeSteps);
        S_arr[i] = path[timeSteps];
        Scap += path[timeSteps];
        if (!isBarrierOption || (shouldBarrierBeHit && isBarrierHit(path, timeSteps)) // downin and upin options - Barrier should be hit and it is hit
            || (!shouldBarrierBeHit && !isBarrierHit(path, timeSteps))                // downout and upout options - Barrier should not be hit and it is not hit
        )
        {
            double priceForPath = exp(-_r * _T) * getPathPayoff(path[timeSteps]);
            V_arr[i] = priceForPath;
            Vcap += priceForPath;
        }
    }

    Scap = Scap / numPaths;
    Vcap = Vcap / numPaths;

    double bnum = 0;
    double bden = 0;
    for (std::size_t i = 0; i < numPaths; i++)
    {
        bnum += (S_arr[i] - Scap) * (V_arr[i] - Vcap);
        bden += (S_arr[i] - Scap) * (S_arr[i] - Scap);
    }
    double bcap = bnum / bden;

    double Wcap = 0;
    for (std::size_t i = 0; i < numPaths; i++)
    {
        Wcap += V_arr[i] - bcap * (S_arr[i] - exp(_r * _T) * _S);
    }
    Wcap = Wcap / numPaths;
    price_return.push_back(Wcap);
    return price_return;
}

std::vector<double> MonteCarlo::price_option_antithetic(std::size_t numPaths, std::size_t timeSteps, bool include_greeks)
{
    MonteCarlo::validateParameters();
    bool isBarrierOption = _option.isBarrierOption();
    bool shouldBarrierBeHit = _type == OptionType::downin || _type == OptionType::upin;

    std::vector<double> price_return;

    std::size_t len = numPaths * timeSteps + 1;
    double *randoms = new double[len];
    double *randomsFlipped = new double[len];
    RandomNumberGenerator *rng = getRandomNumberGenerator();
    for (std::size_t i = 0; i < len; i++)
    {
        double rdm = (*rng)();
        randoms[i] = rdm;
        randomsFlipped[i] = -rdm;
    }

    double *path, *pathFlipped;
    double Vcap = 0;
    double VcapFlipped = 0;
    for (std::size_t i = 0; i < numPaths; i++)
    {
        path = generatePathSDE(&randoms[i * timeSteps], timeSteps);
        pathFlipped = generatePathSDE(&randomsFlipped[i * timeSteps], timeSteps);
        if (!isBarrierOption || (shouldBarrierBeHit && isBarrierHit(path, timeSteps)) // downin and upin options - Barrier should be hit and it is hit
            || (!shouldBarrierBeHit && !isBarrierHit(path, timeSteps))                // downout and upout options - Barrier should not be hit and it is not hit
        )
        {
            Vcap += exp(-_r * _T) * getPathPayoff(path[timeSteps]);
            VcapFlipped += exp(-_r * _T) * getPathPayoff(pathFlipped[timeSteps]);
        }
    }

    Vcap = (Vcap + VcapFlipped) / 2 / numPaths;
    price_return.push_back(Vcap);
    return price_return;
}

std::vector<double> MonteCarlo::price_option_momentmatching(std::size_t numPaths, std::size_t timeSteps, bool include_greeks)
{
    MonteCarlo::validateParameters();
    bool isBarrierOption = _option.isBarrierOption();
    bool shouldBarrierBeHit = _type == OptionType::downin || _type == OptionType::upin;

    std::vector<double> price_return;

    std::size_t len = numPaths * timeSteps + 1;
    double *randoms = new double[len];
    RandomNumberGenerator *rng = getRandomNumberGenerator();
    for (std::size_t i = 0; i < len; i++)
    {
        randoms[i] = (*rng)();
    }

    double *path;
    double *S_arr = new double[numPaths];
    double *V_arr = new double[numPaths];
    double Scap = 0;
    double Vcap = 0;

    for (std::size_t i = 0; i < numPaths; i++)
    {
        path = generatePathSDE(&randoms[i * timeSteps], timeSteps);
        S_arr[i] = path[timeSteps];
        Scap += path[timeSteps];
        if (!isBarrierOption || (shouldBarrierBeHit && isBarrierHit(path, timeSteps)) // downin and upin options - Barrier should be hit and it is hit
            || (!shouldBarrierBeHit && !isBarrierHit(path, timeSteps))                // downout and upout options - Barrier should not be hit and it is not hit
        )
        {
            double priceForPath = exp(-_r * _T) * getPathPayoff(path[timeSteps]);
            V_arr[i] = priceForPath;
            Vcap += priceForPath;
        }
    }

    Scap = Scap / numPaths;
    Vcap = Vcap / numPaths;

    double *S_tilde = new double[numPaths];
    double *V_tilde = new double[numPaths];
    double Stilde = 0;
    double Vtilde = 0;
    for (std::size_t i = 0; i < numPaths; i++)
    {
        S_tilde[i] = S_arr[i] * exp(_r * _T) * _S / Scap;
        if (!isBarrierOption || (shouldBarrierBeHit && isBarrierHit(path, timeSteps)) // downin and upin options - Barrier should be hit and it is hit
            || (!shouldBarrierBeHit && !isBarrierHit(path, timeSteps))                // downout and upout options - Barrier should not be hit and it is not hit
        )
        {
            V_tilde[i] = exp(-_r * _T) * getPathPayoff(S_tilde[i]);
        }
    }
    Stilde = std::accumulate(S_tilde, S_tilde + numPaths, Stilde) / numPaths;
    Vtilde = std::accumulate(V_tilde, V_tilde + numPaths, Vtilde) / numPaths;

    price_return.push_back(Vtilde);
    return price_return;
}

std::vector<double> MonteCarlo::price_option_controlvariatemomentmatching(std::size_t numPaths, std::size_t timeSteps, bool include_greeks)
{
    MonteCarlo::validateParameters();
    bool isBarrierOption = _option.isBarrierOption();
    bool shouldBarrierBeHit = _type == OptionType::downin || _type == OptionType::upin;

    std::vector<double> price_return;

    std::size_t len = numPaths * timeSteps + 1;
    double *randoms = new double[len];
    RandomNumberGenerator *rng = getRandomNumberGenerator();
    for (std::size_t i = 0; i < len; i++)
    {
        randoms[i] = (*rng)();
    }

    double *path;
    double *S_arr = new double[numPaths];
    double *V_arr = new double[numPaths];
    double Scap = 0;
    double Vcap = 0;

    for (std::size_t i = 0; i < numPaths; i++)
    {
        path = generatePathSDE(&randoms[i * timeSteps], timeSteps);
        S_arr[i] = path[timeSteps];
        Scap += path[timeSteps];
        if (!isBarrierOption || (shouldBarrierBeHit && isBarrierHit(path, timeSteps)) // downin and upin options - Barrier should be hit and it is hit
            || (!shouldBarrierBeHit && !isBarrierHit(path, timeSteps))                // downout and upout options - Barrier should not be hit and it is not hit
        )
        {
            double priceForPath = exp(-_r * _T) * getPathPayoff(path[timeSteps]);
            V_arr[i] = priceForPath;
            Vcap += priceForPath;
        }
    }

    Scap = Scap / numPaths;
    Vcap = Vcap / numPaths;

    double *S_tilde = new double[numPaths];
    double *V_tilde = new double[numPaths];
    double Stilde = 0;
    double Vtilde = 0;
    for (std::size_t i = 0; i < numPaths; i++)
    {
        S_tilde[i] = S_arr[i] * exp(_r * _T) * _S / Scap;
        if (!isBarrierOption || (shouldBarrierBeHit && isBarrierHit(path, timeSteps)) // downin and upin options - Barrier should be hit and it is hit
            || (!shouldBarrierBeHit && !isBarrierHit(path, timeSteps))                // downout and upout options - Barrier should not be hit and it is not hit
        )
        {
            V_tilde[i] = exp(-_r * _T) * getPathPayoff(S_tilde[i]);
        }
    }
    Stilde = std::accumulate(S_tilde, S_tilde + numPaths, Stilde) / numPaths;
    Vtilde = std::accumulate(V_tilde, V_tilde + numPaths, Vtilde) / numPaths;

    double bnum = 0;
    double bden = 0;
    for (std::size_t i = 0; i < numPaths; i++)
    {
        bnum += (S_tilde[i] - exp(_r * _T) * _S) * (V_tilde[i] - Vtilde);
        bden += (S_tilde[i] - exp(_r * _T) * _S) * (S_tilde[i] - exp(_r * _T) * _S);
    }
    double bcap = bnum / bden;

    double Wcap = 0;
    for (std::size_t i = 0; i < numPaths; i++)
    {
        Wcap += V_tilde[i] - bcap * (S_tilde[i] - exp(_r * _T) * _S);
    }
    Wcap = Wcap / numPaths;
    price_return.push_back(Wcap);
    return price_return;
}

std::vector<double> MonteCarlo::price_option_Heston(std::size_t numPaths, std::size_t timeSteps, bool include_greeks)
{
    MonteCarlo::validateParameters();
    bool isBarrierOption = _option.isBarrierOption();
    bool shouldBarrierBeHit = _type == OptionType::downin || _type == OptionType::upin;

    std::vector<double> price_return;

    std::size_t len = numPaths * timeSteps + 1;
    double *randomSpots = new double[len];
    RandomNumberGenerator *rng = getRandomNumberGenerator();
    for (std::size_t i = 0; i < len; i++)
    {
        randomSpots[i] = (*rng)();
    }

    double *randomVols = new double[len];
    for (std::size_t i = 0; i < len; i++)
    {
        randomVols[i] = (*rng)();
    }

    double *path;
    double Vcap = 0;
    for (std::size_t i = 0; i < numPaths; i += 1)
    {
        path = generatePathHeston(&randomSpots[i * timeSteps], &randomVols[i * timeSteps], timeSteps);
        if (!isBarrierOption || (shouldBarrierBeHit && isBarrierHit(path, timeSteps)) // downin and upin options - Barrier should be hit and it is hit
            || (!shouldBarrierBeHit && !isBarrierHit(path, timeSteps))                // downout and upout options - Barrier should not be hit and it is not hit
        )
        {
            Vcap += exp(-_r * _T) * getPathPayoff(path[timeSteps]);
        }
    }

    Vcap = Vcap / numPaths;
    price_return.push_back(Vcap);
    return price_return;
}

bool MonteCarlo::isBarrierHit(double *path, std::size_t timeSteps)
{
    double B = _add_params[0];
    bool barrierHit = false;
    for (std::size_t j = 0; j <= timeSteps; j++)
    {
        if (_type == OptionType::downin || _type == OptionType::downout)
        {
            if (path[j] <= B)
            {
                barrierHit = true;
                break;
            }
        }
        else if (_type == OptionType::upin || _type == OptionType::upout)
        {
            if (path[j] >= B)
            {
                barrierHit = true;
                break;
            }
        }
    }

    return barrierHit;
}

double MonteCarlo::getPathPayoff(double lastPathValue)
{
    if (_payoff == OptionPayoff::call)
    {
        return std::max(lastPathValue - _K, 0.0);
    }

    if (_payoff == OptionPayoff::put)
    {
        return std::max(_K - lastPathValue, 0.0);
    }

    throw std::invalid_argument("Payoff not defined for option type");
}