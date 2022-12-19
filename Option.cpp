//
//  Option.cpp
//  MonteCarlo
//
//  Created by Nicolas Buchwalder on 28.11.22.
//

#include <iostream>
#include "Option.hpp"

Option::Option(OptionExercise ex, OptionPayoff payoff, OptionType type, double S, double K, double T, double sigma, double r, double q, DivsTuple divs, std::vector<double> add_params)
    : _ex(ex), _payoff(payoff), _type(type), _S(S), _K(K), _T(T), _sigma(sigma), _r(r), _q(q), _divs(divs), _add_params(add_params), _t(0)
{

    check_barrier_validity();
};

ParamsTuple Option::get_params() const
{
    return std::make_tuple(_ex, _payoff, _type, _S, _K, _T, _sigma, _r, _q, _divs, _add_params);
}

void Option::update_params(){
    
    _d1 = (std::log(_S / _K) + (_r - _q + _sigma * _sigma / 2.) * (_T - _t)) / (_sigma * std::sqrt(_T - _t));
    _d2 = _d1 - _sigma * std::sqrt(_T - _t);

    _zd1 = this->z(_d1);
    _zd2 = this->z(_d2);
    _Nd1 = this->phi(_d1);
    _Nd2 = this->phi(_d2);

    _q_disc = std::exp(-_q * (_T - _t));
    _r_disc = std::exp(-_r * (_T - _t));
}

void Option::set_to_time(double t)
{
    _t = t;
    update_params();
};


void Option::update_price(double S)
{
    _S = S;
    update_params();
};

void Option::update_price_and_expiration(double S, double T)
{
    _S = S;
    _T = T;
    update_params();
};


void Option::check_barrier_validity() const
{
    if (_type != OptionType::vanilla && _add_params.size() == 0)
    {
        throw std::invalid_argument("Option type is not vanilla but no barrier defined.");
    }
}

// CDF of std normal
double Option::phi(double t) const
{
    return std::erfc(-t / std::sqrt(2.)) / 2.;
};

// PDF of std normal
double Option::z(double t) const
{
    return std::exp(-t * t / 2.) / std::sqrt(2. * std::numbers::pi);
}

std::vector<double> Option::price_european(bool includeGreeks)
{
    update_params();
    std::vector<double> res;
    switch (_payoff)
    {
    case call:
        res.push_back(_S * _q_disc * _Nd1 - _K * _r_disc * _Nd2); // price
        if (includeGreeks)
        {
            res.push_back(_q_disc * _Nd1);                                                                                                    // delta
            res.push_back(_q_disc / (_S * _sigma * std::sqrt(_T - _t)) * _zd1);                                                               // gamma
            res.push_back(_S * _q_disc * std::sqrt(_T - _t) * _zd1);                                                                          // vega
            res.push_back(-(_S * _sigma * _q_disc) / (2. * std::sqrt(_T - _t)) * _zd1 + _q * _S * _q_disc * _Nd1 - _r * _K * _r_disc * _Nd2); // theta
            res.push_back(_K * (_T - _t) * _r_disc * _Nd2);                                                                                   // rho
        }
        break;

    case put:
        res.push_back(-_S * _q_disc * (1. - _Nd1) + _K * _r_disc * (1. - _Nd2)); // price
        if (includeGreeks)
        {
            res.push_back(-_q_disc * (1 - _Nd1));                                                                                                         // delta
            res.push_back(_q_disc / (_S * _sigma * std::sqrt(_T - _t)) * _zd1);                                                                           // gamma
            res.push_back(_S * _q_disc * std::sqrt(_T - _t) * _zd1);                                                                                      // vega
            res.push_back(-(_S * _sigma * _q_disc) / (2. * std::sqrt(_T - _t)) * _zd1 - _q * _S * _q_disc * (1 - _Nd1) + _r * _K * _r_disc * (1 - _Nd2)); // theta
            res.push_back(-_K * (_T - _t) * _r_disc * (1 - _Nd2));                                                                                        // rho
        }
        break;
    default:
        break;
    }
    return res;
};

double Option::price_european_withS(double S){
    double old_S = _S;
    update_price(S);
    double new_price = price_european(false)[0];
    update_price(old_S);
    return new_price;
}

double Option::price_european_withSandT(double S, double T){
    double old_S = _S;
    double old_T = _T;
    update_price_and_expiration(S, T);
    double new_price = price_european(false)[0];
    update_price_and_expiration(old_S, old_T);
    return new_price;
}



double Option::calculate_iv(double actualprice, double tol)
{
    double x0 = 0.25;
    double x_old = 0;
    double x_new = x0;
    while (std::abs(x_new - x_old) > tol)
    {
        x_old = x_new;
        _sigma = x_new;
        x_new = x_new - (price_european()[0] - actualprice) / price_european(true)[3];
    }
    return x_new;
}


bool Option::isBarrierOption() const
{
    return _type != OptionType::vanilla;
}

