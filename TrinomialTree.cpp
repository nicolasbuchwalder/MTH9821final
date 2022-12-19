//
//  TrinomialTree.cpp
//  FD
//
//  Created by Nicolas Buchwalder on 19.12.22.
//

#include "TrinomialTree.hpp"

void TrinomialTree::set_payoff(){

    if (_type == vanilla && _payoff == OptionPayoff::call){
        _payoff_now = [=](double S){ return std::max(S - _K, 0.); };
        _price_BS = [=](double S, double T) { return _option.price_european_withSandT(S, T);};
    }
    if (_type == vanilla && _payoff == OptionPayoff::put){
        _payoff_now = [=](double S){ return std::max(_K - S, 0.); };
        _price_BS = [=](double S, double T) { return _option.price_european_withSandT(S, T);};
    }
}

void TrinomialTree::set_terminal_values_standard(int N){
    curr_Vs.resize(2 * N + 1);
    curr_Vs[N] = _payoff_now(_S);
    double S_plus = _S;
    double S_minus = _S;
    for (int i = 1; i <= N; i++){
        S_plus *= _u;
        S_minus *= _d;
        curr_Vs[N+i] = _payoff_now(S_plus);
        curr_Vs[N-i] = _payoff_now(S_minus);
    }
}

void TrinomialTree::set_terminal_values_black_scholes(int N){
    if (_ex == euro){
        set_terminal_values_black_scholes_european(N);
    }
    else {
        set_terminal_values_black_scholes_american(N);
    }
}

void TrinomialTree::set_terminal_values_black_scholes_european(int N){
    curr_Vs.resize(2 * N + 1);
    curr_Vs[N] = _price_BS(_S, _dt);
    double S_plus = _S;
    double S_minus = _S;
    for (int i = 1; i <= N; i++){
        S_plus *= _u;
        S_minus *= _d;
        curr_Vs[N+i] = _price_BS(S_plus, _dt);
        curr_Vs[N-i] = _price_BS(S_minus, _dt);
    }
}

void TrinomialTree::set_terminal_values_black_scholes_american(int N){
    curr_Vs.resize(2 * N + 1);
    curr_Vs[N] = std::max(_price_BS(_S, _dt), _payoff_now(_S));
    double S_plus = _S;
    double S_minus = _S;
    for (int i = 1; i <= N; i++){
        S_plus *= _u;
        S_minus *= _d;
        curr_Vs[N+i] = std::max(_price_BS(S_plus, _dt), _payoff_now(S_plus));
        curr_Vs[N-i] = std::max(_price_BS(S_minus, _dt), _payoff_now(S_minus));
    }
}

void TrinomialTree::backtrack(int N){
    if (_ex == euro){
        backtrack_european(N);
    }
    else {
        backtrack_american(N);
    }
}

void TrinomialTree::backtrack_european(int N){
    for (int j = N - 1; j >= 0; j--){
        if (j == 0){
            past_past_Vs = std::move(past_Vs);
        }
        past_Vs = std::move(curr_Vs);
        for (int i = 0; i <= 2 * j; i++){
            curr_Vs.push_back(_p_up_disc * past_Vs[i] + _p_m_disc * past_Vs[i+1] + _p_down_disc * past_Vs[i+2]);
        }
    }
}

void TrinomialTree::backtrack_american(int N){
    double S_temp;
    for (int j = N - 1; j >= 0; j--){
        if (j == 0){
            past_past_Vs = std::move(past_Vs);
        }
        past_Vs = std::move(curr_Vs);
        S_temp = _S * std::pow(_u, j);
        for (int i = 0; i <= 2 * j; i++){
            curr_Vs.push_back(std::max(_p_up_disc * past_Vs[i] + _p_m_disc * past_Vs[i+1] + _p_down_disc * past_Vs[i+2], _payoff_now(S_temp)));
            S_temp *= _d;
        }
    }
}


std::vector<double> TrinomialTree::price_standard(int N){
    
    update_tree_params(N);
    set_terminal_values_standard(N);
    backtrack(N);
    
    std::vector<double> res;
    res.push_back(curr_Vs[0]);
    
    std::vector<double> gs = greeks();
    res.insert(res.end(), gs.begin(), gs.end());
    
    return res;
}


std::vector<double> TrinomialTree::price_average(int N){
    
    update_tree_params(N + 1);
    set_terminal_values_standard(N + 1);
    backtrack(N + 1);
    double firstprice = curr_Vs[0];
    std::vector<double> gs1 = greeks();
    
    update_tree_params(N);
    set_terminal_values_standard(N);
    backtrack(N);
    double secondprice = curr_Vs[0];
    std::vector<double> gs2 = greeks();
    
    std::vector<double> res;
    res.push_back((firstprice + secondprice) / 2.);
    for (int i = 0; i < gs1.size(); i++){
        res.push_back((gs1[i] + gs2[i]) / 2. );
    }
    return res;
}


std::vector<double> TrinomialTree::price_bbs(int N){
    
    update_tree_params(N);
    set_terminal_values_black_scholes(N);
    backtrack(N - 1);
    
    std::vector<double> res;
    res.push_back(curr_Vs[0]);
    
    std::vector<double> gs = greeks();
    res.insert(res.end(), gs.begin(), gs.end());
    
    clear_contents();
    
    return res;
}


std::vector<double> TrinomialTree::price_bbsr(int N){
    
    update_tree_params(N);
    set_terminal_values_black_scholes(N);
    backtrack(N - 1);
    double firstprice = curr_Vs[0];
    std::vector<double> gs1 = greeks();
    
    update_tree_params(N / 2);
    set_terminal_values_black_scholes(N / 2);
    backtrack(N / 2 - 1);
    double secondprice = curr_Vs[0];
    std::vector<double> gs2 = greeks();
    clear_contents();
    
    std::vector<double> res;
    res.push_back(2 * firstprice - secondprice);
    for (int i = 0; i < gs1.size(); i++){
        res.push_back(2 * gs1[i] - gs2[i]);
    }
    return res;
}


std::vector<double> TrinomialTree::greeks(){
    std::vector<double> past_Ss = {_S * _u, _S * _d};
    std::vector<double> past_past_Ss = {_S * _u * _u, _S, _S * _d * _d};
    double delta = (past_Vs[0] - past_Vs[1]) / (past_Ss[0] - past_Ss[1]);
    double delta_up = (past_past_Vs[0] - past_past_Vs[1]) / (past_past_Ss[0] - past_past_Ss[1]);
    double delta_down = (past_past_Vs[1] - past_past_Vs[2]) / (past_past_Ss[1] - past_past_Ss[2]);
    double gamma = (delta_up - delta_down) / ((past_past_Ss[0] - past_past_Ss[2]) / 2.);
    double theta = (past_past_Vs[1] - curr_Vs[0]) / (2 * _dt);
    return std::vector<double>{delta, gamma, theta};
}

void TrinomialTree::clear_contents(){
    curr_Vs.clear();
}

/*
 PRINT FUNCTIONS
 */


template < typename T >
void TrinomialTree::print(const std::vector<T>& vec) const {
    for (auto elem : vec) {
        std::cout << elem << "\t";
    }
    std::cout << std::endl;
}
template < typename T >
void TrinomialTree::print(const std::vector<std::vector<T>>& mat) const {
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

void TrinomialTree::set_params(Option opt, int init_N, double tol){
    
    _option = opt;
    std::tie(_ex, _payoff, _type, _S, _K, _T, _sigma, _r, _q, _divs, _add_params) = opt.get_params();
    
    if (!_divs.empty()) _q = 0;
    _num_divs = _divs.size();
    
    _init_N = init_N;
    _tol = tol;
    
    update_tree_params(init_N);
    
    set_payoff();
    
}

void TrinomialTree::update_tree_params(int N){
    
    _N = N;
    _dt = _T / _N;
    _u = std::exp(_sigma * std::sqrt(3 * _dt));
    _d = 1. / _u;
    _del = (_r - _q  - _sigma * _sigma / 2.) * std::sqrt(_dt / (12. * _sigma * _sigma));
    _r_disc = std::exp(-_r * _dt);
    _p_up_disc = (1. / 6.  + _del) * _r_disc;
    _p_m_disc = 2. / 3. * _r_disc;
    _p_down_disc = (1. / 6.  - _del) * _r_disc;
    
    clear_contents();
}


std::vector<double> TrinomialTree::price_option_no_tol(const Method& method){
    std::vector<double> res;
    switch (method) {
        case standard:
            res = price_standard(_init_N);
            break;
        case average:
            res = price_average(_init_N);
            break;
        case bbs:
            res = price_bbs(_init_N);
            break;
        case bbsr:
            res = price_bbsr(_init_N);
            break;
        default:
            std::cout << "method not recognised" << std::endl;
            break;
    }
    
    return res;
    
    
}

std::vector<double> TrinomialTree::price_option(const Method& method){
    std::vector<double> res;
    switch (method) {
        case standard: {
            res = price_standard(_init_N);
            int current_N = _init_N * 2;
            double old_price = 0, new_price = res[0];
            while (std::abs(new_price - old_price) < _tol){
                old_price = new_price;
                res = price_standard(current_N);
                new_price = res[0];
                current_N *= 2;
            }
            break;
        }
        case average: {
            res = price_average(_init_N);
            int current_N = _init_N * 2;
            double old_price = 0, new_price = res[0];
            while (std::abs(new_price - old_price) < _tol){
                old_price = new_price;
                res = price_standard(current_N);
                new_price = res[0];
                current_N *= 2;
            }
            break;
        }
        case bbs:{
            res = price_bbs(_init_N);
            int current_N = _init_N * 2;
            double old_price = 0, new_price = res[0];
            while (std::abs(new_price - old_price) < _tol){
                old_price = new_price;
                res = price_bbs(current_N);
                new_price = res[0];
                current_N *= 2;
            }
            break;
        }
        case bbsr: {
            res = price_bbsr(_init_N);
            int current_N = _init_N * 2;
            double old_price = 0, new_price = res[0];
            while (std::abs(new_price - old_price) < _tol){
                old_price = new_price;
                res = price_bbsr(current_N);
                new_price = res[0];
                current_N *= 2;
            }
            break;
        }
        default:
            std::cout << "method not recognised" << std::endl;
            break;
    }
    
    return res;
}
