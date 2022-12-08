//
//  Option.hpp
//  MonteCarlo
//
//  Created by Nicolas Buchwalder on 28.11.22.
//

#ifndef Option_hpp
#define Option_hpp

#include <cmath>
#include <numbers>
#include <tuple>
#include <vector>

enum OptionExercise {
    euro,
    amer
};
enum OptionPayoff {
    call,
    put,
    custom
};
enum OptionType {
    vanilla,
    upin,
    downin,
    upout,
    downout
};

using DivsTuple = std::vector<std::tuple<bool, double, double>>;
using ParamsTuple = std::tuple<OptionExercise, OptionPayoff, OptionType, double, double, double, double, double, double, DivsTuple, std::vector<double>>;

class Option {
   private:
    // option parameters
    OptionPayoff _payoff;  // payoff type
    OptionExercise _ex;    // exercise type
    OptionType _type;      // type of option
    double _S;             // Spot price
    double _K;             // Strike price
    double _T;             // Maturity
    double _sigma;         // Volatility
    double _r;             // Const interest rate
    double _q;             // Compound dividend rate

    DivsTuple _divs;  // dividend
    std::vector<double> _add_params;

    double _t;

    // intermediate values
    double _q_disc;
    double _r_disc;
    double _d1;
    double _d2;
    double _zd1;
    double _zd2;
    double _Nd1;
    double _Nd2;

    // helper function
    double z(double t) const;
    double phi(double t) const;

    void update_vals(double t);
    
    
    
public:
    Option() = default;
    
    Option(OptionExercise ex, OptionPayoff payoff, OptionType type, double S, double K, double T, double sigma, double r, double q, DivsTuple divs, std::vector<double> add_params);

    ParamsTuple get_params() const;

    std::vector<double> price_european(bool includeGreeks) const;
};
#endif /* Option_hpp */
