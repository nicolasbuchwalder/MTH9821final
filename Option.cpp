//
//  Option.cpp
//  MonteCarlo
//
//  Created by Nicolas Buchwalder on 28.11.22.
//

#include "Option.hpp"

Option::Option(OptionExercise ex, OptionPayoff payoff, OptionType type, double S, double K, double T, double sigma, double r, double q, DivsTuple divs, std::vector<double> add_params)
    : _ex(ex), _payoff(payoff), _type(type), _S(S), _K(K), _T(T), _sigma(sigma), _r(r), _q(q), _divs(divs), _add_params(add_params), _t(0)
{
    update_vals(0);
};

ParamsTuple Option::get_params() const
{
    return std::make_tuple(_ex, _payoff, _type, _S, _K, _T, _sigma, _r, _q, _divs, _add_params);
}

void Option::update_vals(double t)
{
    _t = t;

    _d1 = (std::log(_S / _K) + (_r - _q + _sigma * _sigma / 2.) * (_T - _t)) / (_sigma * std::sqrt(_T - _t));
    _d2 = _d1 - _sigma * std::sqrt(_T - _t);

    _zd1 = this->z(_d1);
    _zd2 = this->z(_d2);
    _Nd1 = this->phi(_d1);
    _Nd2 = this->phi(_d2);

    _q_disc = std::exp(-_q * (_T - _t));
    _r_disc = std::exp(-_r * (_T - _t));
};

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

std::vector<double> Option::price_european(bool includeGreeks /*= false*/) const
{
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
