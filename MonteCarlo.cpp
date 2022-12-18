//
//  FiniteDifference.cpp
//  FD
//
//  Created by Nicolas Buchwalder on 07.12.22.
//

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

/**
 * Generates one MC path.
 */
double *MonteCarlo::generatePath(double *randoms, std::size_t timeSteps)
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

std::vector<double> MonteCarlo::price_option(std::size_t numPaths, std::size_t timeSteps, bool include_greeks)
{
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
        path = generatePath(&randoms[i * timeSteps], timeSteps);
        bool barrierHit = false;
        for (std::size_t j = 0; j <= timeSteps; j++)
        {
            if (path[j] <= _add_params[0])
            {
                barrierHit = true;
                break;
            }
        }
        if (!barrierHit)
        {
            Vcap += exp(-_r * _T) * std::max(path[timeSteps] - _K, 0.0);
        }
    }

    price_return.push_back(Vcap / numPaths);
    return price_return;
}