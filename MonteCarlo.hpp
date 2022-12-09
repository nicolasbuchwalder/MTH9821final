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

enum Scheme {
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

class MonteCarlo {
   private:
    OptionPayoff _payoff;  // payoff type
    OptionExercise _ex;    // exercise type
    OptionType _type;      // type of option
    double _S;             // Spot price
    double _K;             // Strike price
    double _T;             // Maturity
    double _sigma;         // Volatility
    double _r;             // Const interest rate
    double _q;             // Compound dividend rate

    std::size_t _num_divs;  // number of discrete dividends
    DivsTuple _divs;        // discrete dividends

    // Extra params
    double _B;

    // list of dividends
    std::vector<double> _tau_divs;
    std::vector<double> _q_divs;

    // computing the greeks
    std::vector<double> greeks();

   public:
    MonteCarlo() = default;
    MonteCarlo(Option opt, std::size_t numPaths, std::size_t timeSteps, RandomNumberGenerator& rng);

    // MC methods
    std::vector<double> price_MC(long long numPaths, long long timeSteps, RandomNumberGenerator rng, bool include_greeks);

    // set all parameters
    void set_params(Option opt, std::size_t M_1, double alpha_1_temp);

    // price option
    std::vector<double> price_option(const Scheme& scheme, bool show_domain, bool include_greeks);

    // print the approximations
    void show_grid();
};

#endif /* MonteCarlo_hpp */
