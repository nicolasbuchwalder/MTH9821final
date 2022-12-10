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

enum MCSchemes
{
    standardMC,
    hestonMC
};

/*
enum RandomGeneratorMethod {
    LinearCongruential,
    InverseTransform,
    AcceptanceRejection,
    BoxMuller
}
*/

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

    // list of dividends
    std::vector<double> _tau_divs;
    std::vector<double> _q_divs;

    // computing the greeks
    std::vector<double> greeks();

public:
    MonteCarlo() = default;
    MonteCarlo(Option opt, std::size_t numPaths, std::size_t timeSteps, RandomNumberGenerator &rng);

    // MC methods
    std::vector<double> price_option(std::size_t numPaths, std::size_t timeSteps, std::vector<double> randoms, bool include_greeks);
    std::vector<double> price_MC_main(std::vector<std::size_t> numPathsVec, std::vector<std::size_t> timeStepsVec, RandomNumberGenerator &random, bool include_greeks);
};

#endif /* MonteCarlo_hpp */
