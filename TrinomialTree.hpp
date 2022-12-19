//
//  TrinomialTree.hpp
//  FD
//
//  Created by Nicolas Buchwalder on 19.12.22.
//

#ifndef TrinomialTree_hpp
#define TrinomialTree_hpp

#include "Option.hpp"
#include "BinomialTree.hpp"

#include <iostream>


class TrinomialTree {

private:
    
    Option _option;
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
    std::size_t _num_divs;
    std::vector<double> _add_params;
   
    // tree params
    int _init_N;
    double _tol;
    
    int _N;
    double _dt;
    double _u, _d;
    double _del;
    double _r_disc;
    double _p_up_disc, _p_m_disc, _p_down_disc;
    
    std::function<double (double)> _payoff_now;
    std::function<double (double, double)> _price_BS;
    
    std::vector<double> curr_Vs;
    std::vector<double> past_Vs;
    std::vector<double> past_past_Vs;
    
    void update_tree_params(int N);
    void clear_contents();
    
    void set_payoff();
    
    void set_terminal_values_standard(int N);
    void set_terminal_values_black_scholes(int N);
    void set_terminal_values_black_scholes_european(int N);
    void set_terminal_values_black_scholes_american(int N);
    
    void backtrack(int N);
    void backtrack_european(int N);
    void backtrack_american(int N);
    
    
    std::vector<double> price_standard(int N);
    std::vector<double> price_average(int N);
    std::vector<double> price_bbs(int N);
    std::vector<double> price_bbsr(int N);
    
    std::vector<double> greeks();
    

public:
    
    TrinomialTree() = default;
    
    // set all parameters
    void set_params(Option opt, int init_N, double tol);
    
    std::vector<double> price_option_no_tol(const Method& method);
    std::vector<double> price_option(const Method& method);
    
    
    // print functions
    template < typename T >
    void print(const std::vector<T>& vec) const;
    template < typename T >
    void print(const std::vector<std::vector<T>>& mat) const;
    
};

#endif /* TrinomialTree_hpp */
