//
//  FiniteDifference.cpp
//  FD
//
//  Created by Nicolas Buchwalder on 07.12.22.
//

#include "FiniteDifference.hpp"

#include "LinearSolver.hpp"
#include "Decomposer.hpp"
#include "IterativeSolver.hpp"


/*
 DOMAIN BUILDERS
 */

void FiniteDifference::set_domain(){
    
    _tau_final_ = _T * _sigma * _sigma / 2.;
    
    for (auto t: _divs){
        _tau_divs.push_back((_T - std::get<1>(t)) * _sigma * _sigma / 2.);
        _q_divs.push_back(std::get<2>(t));
    }
    
    _x_l = std::log(_S / _K) + (_r - _sigma * _sigma / 2.) * _T - 3. * _sigma * std::sqrt(_T);
    _x_r = std::log(_S / _K) + (_r - _sigma * _sigma / 2.) * _T + 3. * _sigma * std::sqrt(_T);
    
    if (_type == OptionType::downout){
        _x_l = std::log(_add_params[0] / _K);
    }
    
    if (_type == OptionType::upout){
        _x_r = std::log(_add_params[0] / _K);
    }
    
}

void FiniteDifference::set_discretisation(){
    double x_compute = std::log(_S / _K);
    for (auto q : _q_divs){x_compute += (1 - q);};
    
    if 
    
    for (int i = 1; i < _num_divs; i++){
        _dtaus.push_back(_tau_divs[i] / _Ms[i]);
        _dxs.push_back(std::sqrt(_dtaus[i] / _alphas[i]));
        
        
    }
    double N_left = std::ceil((x_bar_compute - x_l_) / dx);
    double N_right = std::ceil((x_r_ - x_bar_compute) / dx);
    double N = N_left + N_right;
    
    //std::cout << "before:" << x_l_ << ", " << x_r_ << std::endl;
    x_l_ = x_bar_compute - N_left * dx;
    x_r_ = x_bar_compute + N_right * dx;
    //std::cout << "after:" << x_l_ << ", " << x_r_ << std::endl;
    
    x_l_new_ = x_l_ - std::log(1 - q_);
    x_r_new_ = x_r_ - std::log(1 - q_);
    
    double dtau_2 = alpha_temp_ * dx * dx;
    double M_2 = std::ceil((tau_final_ - tau_div_) / dtau_2);
    dtau_2 = (tau_final_ - tau_div_) / M_2;
    double alpha_2 = dtau_2 / (dx * dx);
    
}

void FiniteDifference::build_mesh(){
    // building the mesh on x
    _x_mesh.push_back(_x_l);
    for (std::size_t i = 0; i < _N; i++) {
        _x_mesh.push_back(_x_mesh.back() + _dxs.front());
    };
    
    // generating first layer of approximations
    _u_mesh.push_back(std::vector<double>(_N));
    std::transform(_x_mesh.cbegin(), _x_mesh.cend(), _u_mesh.back().begin(), _boundary_tau_0);
    
}


//

void FiniteDifference::set_boundaries(){
    // VANILLA EUROPEAN CALL
    if (_payoff == OptionPayoff::call && _type == vanilla){
        _boundary_tau_0 = [=](double x)->double {return std::max(_K * std::exp(_a * x) * (std::exp(x) - 1.), 0.);};
        _boundary_x_l = [](double tau, double x)->double {return 0.;};
        _boundary_x_r = [=](double tau, double x)->double { return _K * std::exp(_a * x + _b * tau) * (std::exp(_x_r) - std::exp(-2. * _r * tau / (_sigma * _sigma)));};
    }
    
    
    
}


/*
 STEP ADVANCE
 */

void FiniteDifference::advance_expl(double alpha, double tau, double xl, double xr){
    
    std::vector<double> new_u_mesh;
    
    // left boundary
    new_u_mesh.push_back(_boundary_x_l(tau, xl));
    
    // middle values
    for (std::size_t pos = 1; pos < _N - 1; pos++) {
        new_u_mesh.push_back(alpha * _u_mesh.back()[pos - 1] + (1. - 2. * alpha) * _u_mesh.back()[pos] + alpha * _u_mesh.back()[pos + 1]);
    }
    
    // right boundary
    new_u_mesh.push_back(_boundary_x_r(tau, xr));
    
    _u_mesh.push_back(new_u_mesh);
}


void FiniteDifference::advance_impl(double alpha, double tau, double xl, double xr, const mat& L, const mat& U){
    
    std::vector<double> new_u_mesh(_u_mesh.back().size());
    
    vec b(_u_mesh.back().size() - 2);
    std::copy(_u_mesh.back().cbegin() + 1, _u_mesh.back().cend() - 1, b.begin());
    
    // add boundary conditions
    b(0) += _boundary_x_l(tau, xl) * alpha;
    b(b.size() - 1) += _boundary_x_r(tau, xr) * alpha;
    
    // solve linear system
    LinearSolver substituter;
    b = substituter.ForwardSub(L, b);
    b = substituter.BackwardSub(U, b);
    
    // assign to new u_mesh
    std::copy(b.cbegin(), b.cend(), new_u_mesh.begin() + 1);
    
    // add boundary condition
    *(new_u_mesh.begin()) = _boundary_x_l(tau, xl);
    *(new_u_mesh.rbegin()) = _boundary_x_r(tau, xr);
    
    _u_mesh.push_back(new_u_mesh);
}








/*
 SUB DOMAIN PRICERS
 */

void FiniteDifference::compute_sub_domain_expl(){
    
}






/*
 GLOBAL PRICERS
 */

std::vector<double> FiniteDifference::price_expl(bool show_domain, bool include_greeks){
    
    std::vector<double> res;
    
    
    
    return res;
}






/*
 PUBLIC FUNCTIONS
 */

void FiniteDifference::set_params(Option opt, std::size_t M_1, double alpha_1_temp){
    
    std::tie(_ex, _payoff, _type, _S, _K, _T, _sigma, _r, _q, _divs, _add_params) = opt.get_params();
    
    // making sure no continuous compounding if discrete dividends
    if (!_divs.empty()) _q = 0;
    _num_divs = _divs.size();
    
    
    auto num_sub_domains = _divs.size() + 1;
    
    _Ms.assign(num_sub_domains, 0.);
    _Ms[0] = M_1;
    _alphas.assign(num_sub_domains, 0.);
    _alphas[0] = alpha_1_temp;
    
    
    // heat coefficients
    _a = (_r - _q) / (_sigma * _sigma) - .5;
    _b = (_a + 1.) * (_a + 1.) + 2 * _q / (_sigma * _sigma);
    
}

std::vector<double> FiniteDifference::price_option(const Scheme& scheme, bool show_domain, bool include_greeks){
    
    std::vector<double> res;
    switch (scheme){
        case eul_expl:
            res = price_expl(include_greeks);
            break;
        default:
            break;
    };
    return res;
};


void FiniteDifference::show_grid(){
    for (auto& row : _u_mesh){
        for (auto& elem : row){
            std::cout << elem << "\t";
        };
        std::cout << std::endl;
    }
}







