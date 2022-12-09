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

MonteCarlo::MonteCarlo(Option option, std::size_t numPaths, std::size_t timeSteps, RandomNumberGenerator &rng)
{
    std::tie(_ex, _payoff, _type, _S, _K, _T, _sigma, _r, _q, _divs, _add_params) = option.get_params();
    _option = option;
}

std::vector<double> MonteCarlo::price_option(std::size_t numPaths, std::size_t timeSteps, std::vector<double> randoms, bool include_greeks)
{
    double dt = _T / timeSteps;
    for (std::size_t i = 0; i < numPaths; i++)
    {
    }
}

std::vector<double> MonteCarlo::price_MC_main(std::vector<std::size_t> numPathsVec, std::vector<std::size_t> timeStepsVec, RandomNumberGenerator &random, bool include_greeks)
{
    std::size_t N = numPathsVec.back() * timeStepsVec.back();
    std::vector<double> randoms;
    for (std::size_t i = 0; i < N; i++)
    {
        randoms.push_back(random());
    }

    std::size_t len = numPathsVec.size();

    for (int i = 0; i < len; i++)
    {
        std::vector<double> price_BS = _option.price_european();

        std::size_t numPaths = numPathsVec[i];
        std::size_t timeSteps = timeStepsVec[i];
        std::vector<double> price_res = price_option(numPaths, timeSteps, randoms, include_greeks);

        std::size_t timeSteps2 = ceil(pow(N, 1.0 / 3) * pow(_T, 2.0 / 3));
        std::size_t numPaths2 = floor((double)N / timeSteps);
        std::vector<double> price_res2 = price_option(numPaths2, timeSteps2, randoms, include_greeks);

        std::cout << N << " " << timeSteps << " " << numPaths << " " << price_res[0] << abs(price_BS[0] - price_res[0]) << " " << timeSteps2 << " " << numPaths2 << " " << price_res2[0] << " " << abs(price_BS[0] - price_res2[0]) << std::endl;
    }
}
