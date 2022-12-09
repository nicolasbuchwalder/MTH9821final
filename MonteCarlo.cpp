//
//  FiniteDifference.cpp
//  FD
//
//  Created by Nicolas Buchwalder on 07.12.22.
//

#include "MonteCarlo.hpp"

#include "Decomposer.hpp"
#include "FiniteDifference.hpp"
#include "IterativeSolver.hpp"
#include "LinearSolver.hpp"
#include "RandomNumberGenerator.hpp"

MonteCarlo::MonteCarlo(Option opt, std::size_t numPaths, std::size_t timeSteps, RandomNumberGenerator& rng) {
    std::tie(_ex, _payoff, _type, _S, _K, _T, _sigma, _r, _q, _divs, _add_params) = opt.get_params();

    OptionExercise _ex;  // exercise type
    OptionType _type;    // type of option
    double _S;           // Spot price
    double _K;           // Strike price
    double _T;           // Maturity
    double _sigma;       // Volatility
    double _r;           // Const interest rate
    double _q;           // Compound dividend rate

    std::size_t _num_divs;  // number of discrete dividends
    DivsTuple _divs;        // discrete dividends

    // Extra params
    double _B;
}
std::vector<double> price_MC(long long numPaths, long long timeSteps, RandomNumberGenerator& random, bool include_greeks) {
    long long N = numPaths * timeSteps;
    double* randoms = new double[N];
    for (long long i = 0; i < N; i++) {
        randoms[i] = random();
    }
}
