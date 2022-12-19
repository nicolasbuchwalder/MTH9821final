//
//  BinomialTree.cpp
//  FD
//
//  Created by Nicolas Buchwalder on 09.12.22.
//

#include "BinomialTree.hpp"

void BinomialTree::set_payoff(){
    
    if (_type == vanilla && _payoff == OptionPayoff::call && _ex == OptionExercise::euro){
        _payoff_mat = [=](double S){ return };
        _early_exercise = [=](double S, double t){ return -std::numeric_limits<double>::max(); };
        
    }
    if (_type == vanilla && _payoff == OptionPayoff::call && _ex == OptionExercise::euro){
        _payoff_function = [=](double S, double t){ return -std::numeric_limits<double>::max(); };
    }
}

void BinomialTree::backtrack(){
    
}

/*
 PRINT FUNCTIONS
 */


template < typename T >
void BinomialTree::print(const std::vector<T>& vec) const {
    for (auto elem : vec) {
        std::cout << elem << "\t";
    }
    std::cout << std::endl;
}
template < typename T >
void BinomialTree::print(const std::vector<std::vector<T>>& mat) const {
    for (auto vec : mat) {
        for (auto elem : vec){
            std::cout << elem << "\t";
        }
        std::cout << std::endl;
    }
}






/*
 PUBLIC METHODS
 */

void BinomialTree::set_params(Option opt, std::size_t N){
    
    _option = opt;
    std::tie(_ex, _payoff, _type, _S, _K, _T, _sigma, _r, _q, _divs, _add_params) = opt.get_params();
    
    if (!_divs.empty()) _q = 0;
    _num_divs = _divs.size();
    
    _N = N;
    _dt = _T / _N;
    _u = std::exp(_sigma * std::sqrt(_dt));
    _d = 1. / _u;
    _p = (std::exp((_r - _q) * _dt) - _d) / (_u - _d);
    _r_disc = std::exp(-_r * _dt);
    _p_disc = _p * _r_disc;
    _q_disc = (1 - _p) * _r_disc;
    
}



std::vector<double> BinomialTree::price_option(const Method& method){
    return std::vector<double>();
}
