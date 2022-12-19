//
//  BinomialTree.hpp
//  FD
//
//  Created by Nicolas Buchwalder on 09.12.22.
//

#ifndef BinomialTree_hpp
#define BinomialTree_hpp

#include "Option.hpp"

#include <iostream>

enum Method {
    standard,
    average,
    bbs,
    bbsr
};

class BinomialTree {

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
    std::size_t _N;
    double _dt;
    double _u, _d;
    double _p;
    double _r_disc;
    double _p_disc, _q_disc;
    
    std::function<double (double)> _payoff_mat;
    std::function<double (double)> _early_exercise;
    
    std::vector<std::vector<double>> V;
    
    void set_payoff();
    
    void set_last_column();
    
    void backtrack();
    
    std::vector<double> greeks();

public:
    
    BinomialTree() = default;
    
    // set all parameters
    void set_params(Option opt, std::size_t N);
    
    
    std::vector<double> price_option(const Method& method);
    
    
    // print functions
    template < typename T >
    void print(const std::vector<T>& vec) const;
    template < typename T >
    void print(const std::vector<std::vector<T>>& mat) const;
    
};
#endif /* BinomialTree_hpp */
