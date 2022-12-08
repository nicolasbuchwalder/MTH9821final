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

void FiniteDifference::build_domain(){
    set_domain();
    set_discretisation();
    set_boundaries();
    build_mesh();
}


void FiniteDifference::set_domain(){
    
    _tau_final = _T * _sigma * _sigma / 2.;
    
    for (auto t: _divs){
        _tau_divs.push_back((_T - std::get<1>(t)) * _sigma * _sigma / 2.);
        _q_divs.push_back(std::get<2>(t));
    }
    _tau_divs.push_back(_tau_final);
    
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
    
    _x_compute = std::log(_S / _K);
    
    // no dividends
    if (_num_divs == 0){
        
        // down out
        if (_type == OptionType::downout){
            _dtaus.push_back(_tau_final / _Ms.front());
            _dxs.push_back(std::sqrt(_dtaus.front() / _alpha_temp));
            std::size_t N_left = std::floor((_x_compute - _x_l) / _dxs.front());
            _x_compute_idx = std::make_pair(N_left, N_left);
            _dxs.front() = (_x_compute - _x_l) / N_left;
            _alphas.push_back(_dtaus.front() / (_dxs.front() * _dxs.front()));
            std::size_t N_right = std::ceil((_x_r - _x_compute) / _dxs.front());
            _N = N_left + N_right;
        }
        
        // up out
        else if (_type == OptionType::upout){
            _dtaus.push_back(_tau_final / _Ms.front());
            _dxs.push_back(std::sqrt(_dtaus.front() / _alpha_temp));
            std::size_t N_right = std::floor((_x_r - _x_compute) / _dxs.front());
            _dxs.front() = (_x_r - _x_compute) / N_right;
            _alphas.push_back(_dtaus.front() / (_dxs.front() * _dxs.front()));
            std::size_t N_left = std::ceil((_x_compute - _x_l) / _dxs.front());
            _x_compute_idx = std::make_pair(N_left, N_left);
            _N = N_left + N_right;
            
        }
        // plain vanilla and barrier in
        else {
            _dtaus.push_back(_tau_final / _Ms.front());
            _N = std::floor((_x_r - _x_l) / std::sqrt(_dtaus.front() / _alpha_temp));
            _dxs.push_back((_x_r - _x_l) / _N);
            _alphas.push_back(_dtaus.front() / (_dxs.front() * _dxs.front()));
            // get index in between value to compute
            double x = _x_l; std::size_t idx = 0;
            while (x < _x_compute){
                x += _dxs.front();
                idx += 1;
            }
            _x_compute_idx = std::make_pair(idx - 1, idx);
            
        };
    }
    // if discrete dividends
    else {
        double x_bar_compute = _x_compute;
        for (auto q : _q_divs){x_bar_compute += (1 - q);};
        double dtau_1 = _tau_divs.front() / _Ms.front();
        _dtaus.push_back(dtau_1);
        double dx = std::sqrt(dtau_1 / _alphas.front());
        _dxs.push_back(dx);
        double N_left = std::ceil((x_bar_compute - _x_l) / dx);
        double N_right = std::ceil((_x_r - x_bar_compute) / dx);
        _x_compute_idx = std::make_pair(N_left, N_left);
        _N = N_left + N_right;
        
        for (int i = 0; i < _num_divs; i++){
            _dtaus.push_back(_alpha_temp * _dxs.front() * _dxs.front());
            _Ms.push_back(std::ceil((_tau_divs[i+1] - _tau_divs[i]) / _dtaus.back()));
            _alphas.push_back(_dtaus.front() / (_dxs.back() * _dxs.back()));
        }
    }
};

void FiniteDifference::build_mesh(){
    // building the mesh on x
    _x_mesh.push_back(_x_l);
    for (std::size_t i = 0; i < _N; i++) {
        _x_mesh.push_back(_x_mesh.back() + _dxs.front());
    };
    // generating first layer of approximations
    _u_mesh.push_back(std::vector<double>(_N+1));
    std::transform(_x_mesh.begin(), _x_mesh.end(), _u_mesh.back().begin(), _boundary_tau_0);
}


//

void FiniteDifference::set_boundaries(){
    // VANILLA EUROPEAN CALL
    if (_payoff == OptionPayoff::call && _type == vanilla){
        _boundary_tau_0 = [=](double x)->double {return std::max(_K * std::exp(_a * x) * (std::exp(x) - 1.), 0.);};
        _boundary_x_l = [](double x, double tau)->double {return 0.;};
        _boundary_x_r = [=](double x, double tau)->double { return _K * std::exp(_a * x + _b * tau) * (std::exp(_x_r) - std::exp(-2. * _r * tau / (_sigma * _sigma)));};
    }
    if (_payoff == OptionPayoff::call && _type == downout){
        _boundary_tau_0 = [=](double x)->double { return std::max(_K * std::exp(_a * x) * (std::exp(x) - 1.), 0.);};
        _boundary_x_l = [](double x, double tau)->double { return 0.;};
        _boundary_x_r = [=](double x, double tau)->double { return _K * std::exp(_a * x + _b * tau) * (std::exp(x - 2. * _q * tau / (_sigma * _sigma) ) - std::exp(-2. * _r * tau / (_sigma * _sigma)));
        };
    }
}


/*
 STEP ADVANCE
 */

void FiniteDifference::advance_expl(double alpha, double tau, double xl, double xr){
    
    std::vector<double> new_u_mesh;
    // left boundary
    new_u_mesh.push_back(_boundary_x_l(xl, tau));
    
    // middle values
    for (std::size_t pos = 1; pos < _N; pos++) {
        new_u_mesh.push_back(alpha * _u_mesh.back()[pos - 1] + (1. - 2. * alpha) * _u_mesh.back()[pos] + alpha * _u_mesh.back()[pos + 1]);
    }
    
    // right boundary
    new_u_mesh.push_back(_boundary_x_r(xr, tau));
    
    _u_mesh.push_back(new_u_mesh);
}


void FiniteDifference::advance_impl(double alpha, double tau, double xl, double xr, const mat& L, const mat& U){
    
    std::vector<double> new_u_mesh(_u_mesh.back().size());
    
    vec b(_u_mesh.back().size() - 2);
    std::copy(_u_mesh.back().cbegin() + 1, _u_mesh.back().cend() - 1, b.begin());
    
    // add boundary conditions
    b(0) += _boundary_x_l(xl, tau) * alpha;
    b(b.size() - 1) += _boundary_x_r(xr, tau) * alpha;
    
    // solve linear system
    LinearSolver substituter;
    b = substituter.ForwardSub(L, b);
    b = substituter.BackwardSub(U, b);
    
    // assign to new u_mesh
    std::copy(b.cbegin(), b.cend(), new_u_mesh.begin() + 1);
    
    // add boundary condition
    *(new_u_mesh.begin()) = _boundary_x_l(xl, tau);
    *(new_u_mesh.rbegin()) = _boundary_x_r(xr, tau);
    
    _u_mesh.push_back(new_u_mesh);
}








/*
 SUB DOMAIN PRICERS
 */

void FiniteDifference::compute_sub_domain_expl(std::size_t sub, double start_tau){
    for (int i = 1; i <= _Ms[sub]; i++){
        advance_expl(_alphas[sub], start_tau + _dtaus[sub] * i, _x_l, _x_r);
    }
    
}






/*
 GLOBAL PRICERS
 */

std::vector<double> FiniteDifference::price_expl(bool include_greeks){
    
    for (int sub = 0; sub <= _num_divs; sub++){
        
        compute_sub_domain_expl(sub, 0);
    }
    
    std::vector<double> res;
    res.push_back(approximate());
    
    return res;
}






/*
 SUPPORT
 */

double FiniteDifference::convert_to_v(double x, double tau, double u) const{
    return std::exp(-_a * x -_b * tau) * u;
}

double FiniteDifference::approximate(){
    // linear interpolation
    double u_approx = ((_x_mesh[_x_compute_idx.second] - _x_compute) * _u_mesh.back()[_x_compute_idx.first] + (_x_compute - _x_mesh[_x_compute_idx.first]) * _u_mesh.back()[_x_compute_idx.second]) / _dxs.back();
    return std::exp(-_a * _x_compute - _b * _tau_final) * u_approx;
    
}


/*
 PRINT FUNCTIONS
 */


template < typename T >
void FiniteDifference::print(const std::vector<T>& vec) const {
    for (auto elem : vec) {
        std::cout << elem << "\t";
    }
    std::cout << std::endl;
}
template < typename T >
void FiniteDifference::print(const std::vector<std::vector<T>>& mat) const {
    for (auto vec : mat) {
        for (auto elem : vec){
            std::cout << elem << "\t";
        }
        std::cout << std::endl;
    }
}



/*
 PUBLIC FUNCTIONS
 */

void FiniteDifference::set_params(Option opt, std::size_t M_1, double alpha_temp){
    
    _option = opt;
    std::tie(_ex, _payoff, _type, _S, _K, _T, _sigma, _r, _q, _divs, _add_params) = opt.get_params();
    
    // making sure no continuous compounding if discrete dividends
    if (!_divs.empty()) _q = 0;
    _num_divs = _divs.size();
    
    _alpha_temp = alpha_temp;
    
    _Ms.clear();
    _Ms.push_back(M_1);
    _alphas.clear();
    
    
    // heat coefficients
    _a = (_r - _q) / (_sigma * _sigma) - .5;
    _b = (_a + 1.) * (_a + 1.) + 2 * _q / (_sigma * _sigma);
    
}

std::vector<double> FiniteDifference::price_option(const Scheme& scheme, bool include_greeks){
    
    build_domain();
    
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

void FiniteDifference::show_domain_params(){
    std::cout << "N: " << _N << std::endl;
    std::cout << "x_l: " << _x_l << std::endl;
    std::cout  <<  "x_r: " << _x_r << std::endl;
    std::cout << "tau_final: " << _tau_final << std::endl;
    std::cout << "tau divs: "; print(_tau_divs);
    std::cout << "dxs: "; print(_dxs);
    std::cout << "Ms: "; print(_Ms);
    std::cout << "dtaus: "; print(_dtaus);
    std::cout << "alphas: "; print(_alphas);
};

void FiniteDifference::show_grid(bool convert){
    
    std::cout << std::endl << "grid:" << std::endl;
     
    if (convert){
        
        double tau = 0;
        double x = _x_l;
        
        for (int i = 0; i < _num_divs + 1; i++ ){
            
            for (int j = 0; j < _Ms[i]; j++){
                
                x = _x_l;
                for (int k = 0; k < _N; k++){
                    std::cout << convert_to_v(x, tau, _u_mesh[j][k]) << "\t";
                    x += _dxs[i];
                }
                tau += _dtaus[i];
                std::cout << std::endl;
            }
        }
    }
    else {
        print(_u_mesh);
    };
    
}







