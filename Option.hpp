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

enum OptionExercise
{
    euro,
    amer
};
enum OptionPayoff
{
    call,
    put,
    custom
};
enum OptionType
{
    vanilla,
    upin,
    downin,
    upout,
    downout
};

// std::tuple<false for fixed and true for proportional, ex-dividend time, value of dividend (can be fixed dollar amount or proportional to stock price)>
using DivsTuple = std::vector<std::tuple<bool, double, double>>;
using ParamsTuple = std::tuple<OptionExercise, OptionPayoff, OptionType, double, double, double, double, double, double, DivsTuple, std::vector<double>>;

class Option
{
private:
    // option parameters
    OptionPayoff _payoff; // payoff type
    OptionExercise _ex;   // exercise type
    OptionType _type;     // type of option
    double _S;            // Spot price
    double _K;            // Strike price
    double _T;            // Maturity
    double _sigma;        // Volatility
    double _r;            // Const interest rate
    double _q;            // Compound dividend rate

    DivsTuple _divs; // dividend
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

    void update_params();
    void check_barrier_validity() const;
    
public:
    Option() = default;

    Option(OptionExercise ex, OptionPayoff payoff, OptionType type, double S, double K, double T, double sigma, double r, double q, DivsTuple divs, std::vector<double> add_params);

    ParamsTuple get_params() const;
    
    void set_to_time(double t);
    void update_price(double S);
    void update_price_and_expiration(double S, double T);
    
    std::vector<double> price_european(bool includeGreeks = false);
    
    double price_european_withS(double S);
    double price_european_withSandT(double S, double T);
    
    
    double calculate_iv(double actualprice, double tol);

    // helper function
    double z(double t) const;
    double phi(double t) const;
    bool isBarrierOption() const;
    
};
#endif /* Option_hpp */
