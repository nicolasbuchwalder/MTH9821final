//
//  MonteCarlo.hpp
//  FD
//
//  Created by Nicolas Buchwalder on 07.12.22.
//

#ifndef MonteCarlo_hpp
#define MonteCarlo_hpp

#include <array>
#include <functional>
#include <iostream>
#include <tuple>
#include <vector>

#include "EigenCommonHeader.h"
#include "Option.hpp"
#include "RandomNumberGenerator.hpp"

class MonteCarlo
{
private:
    Option _option;
    OptionPayoff _payoff; // payoff type
    OptionExercise _ex;   // exercise type
    OptionType _type;     // type of option
    double _S;            // Spot price
    double _K;            // Strike price
    double _T;            // Maturity
    double _sigma;        // Volatility
    double _r;            // Const interest rate
    double _q;            // Compound dividend rate

    std::size_t _num_divs; // number of discrete dividends
    DivsTuple _divs;       // discrete dividends

    // Extra params
    std::vector<double> _add_params;

    // Random Number Generator
    RandomNumberGeneratorEnum _rng;

    // list of dividends
    std::vector<double> _tau_divs;
    std::vector<double> _q_divs;

    // computing the greeks
    std::vector<double> greeks();

public:
    MonteCarlo() = default;
    MonteCarlo(Option opt, RandomNumberGeneratorEnum rng);

    // MC methods
    double *generatePath(double *randoms, std::size_t timeSteps);
    std::vector<double> price_option(std::size_t numPaths, std::size_t timeSteps, bool include_greeks);
    RandomNumberGenerator *getRandomNumberGenerator();
    void validateParameters();
    bool isBarrierHit(double *path, std::size_t timeSteps);
    double getPathPayoff(double s);
};

#endif /* MonteCarlo_hpp */
