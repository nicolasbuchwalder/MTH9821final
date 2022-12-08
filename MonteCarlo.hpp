//
//  MonteCarlo.hpp
//  FD
//
//  Created by Nicolas Buchwalder on 07.12.22.
//

#ifndef MonteCarlo_hpp
#define MonteCarlo_hpp

#include "EigenCommonHeader.h"
#include "Option.hpp"
#include "RandomNumberGenerator.hpp"

#include <functional>
#include <vector>
#include <tuple>
#include <array>
#include <iostream>

enum Scheme {
    standardMC,
    hestonMC
};

class MonteCarlo {

private:
    
    OptionPayoff _payoff;   // payoff type
    OptionExercise _ex;     // exercise type
    OptionType _type;       // type of option
    double _S;              // Spot price
    double _K;              // Strike price
    double _T;              // Maturity
    double _sigma;          // Volatility
    double _r;              // Const interest rate
    double _q;              // Compound dividend rate
    
    std::size_t _num_divs;  // number of discrete dividends
    DivsTuple _divs;        // discrete dividends
    
    // Extra params
    double _B;
    
    // list of dividends
    std::vector<double> _tau_divs;
    std::vector<double> _q_divs;
    
    void advance_expl(double alpha, double tau, double xl, double xr);
    std::vector<double> price_expl(bool show_domain, bool include_greeks);
    
    // computing the greeks
    std::vector<double> greeks();
    
    
public:
    
    MonteCarlo() = default;
    
    // set all parameters
    void set_params(Option opt, std::size_t M_1, double alpha_1_temp);
    
    // price option
    std::vector<double> price_option(const Scheme& scheme, bool show_domain, bool include_greeks);
    
    // print the approximations
    void show_grid();
    
};

#endif /* MonteCarlo_hpp */
