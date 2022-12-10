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
    
    _tau_final = _T * _sigma * _sigma / 2.;
    
    for (auto t: _divs){
        _tau_divs.push_back((_T - std::get<1>(t)) * _sigma * _sigma / 2.);
        _q_divs.push_back(std::get<2>(t));
    }
    _tau_divs.push_back(_tau_final);
    
    _x_ls.push_back(std::log(_S / _K) + (_r - _q - _sigma * _sigma / 2.) * _T - 3. * _sigma * std::sqrt(_T));
    _x_rs.push_back(std::log(_S / _K) + (_r - _q - _sigma * _sigma / 2.) * _T + 3. * _sigma * std::sqrt(_T));
    
    if (_type == OptionType::downout){
        _x_ls[0] = std::log(_add_params[0] / _K);
    }
    
    if (_type == OptionType::upout){
        _x_rs[0] = std::log(_add_params[0] / _K);
    }
    
}

void FiniteDifference::set_discretisation(){
    
    _target_x = std::log(_S / _K);
    
    // no dividends
    if (_num_divs == 0){
        
        // down out
        if (_type == OptionType::downout){
            _dtaus.push_back(_tau_final / _Ms.front());
            _dxs.push_back(std::sqrt(_dtaus.front() / _alpha_temp));
            std::size_t N_left = std::floor((_target_x - _x_ls.front()) / _dxs.front());
            _target_idx = std::make_pair(N_left, N_left);
            _dxs.front() = (_target_x - _x_ls.front()) / N_left;
            _alphas.push_back(_dtaus.front() / (_dxs.front() * _dxs.front()));
            std::size_t N_right = std::ceil((_x_rs.front() - _target_x) / _dxs.front());
            _N = N_left + N_right;
        }
        
        // up out
        else if (_type == OptionType::upout){
            _dtaus.push_back(_tau_final / _Ms.front());
            _dxs.push_back(std::sqrt(_dtaus.front() / _alpha_temp));
            std::size_t N_right = std::floor((_x_rs.front() - _target_x) / _dxs.front());
            _dxs.front() = (_x_rs.front() - _target_x) / N_right;
            _alphas.push_back(_dtaus.front() / (_dxs.front() * _dxs.front()));
            std::size_t N_left = std::ceil((_target_x - _x_ls.front()) / _dxs.front());
            _target_idx = std::make_pair(N_left, N_left);
            _N = N_left + N_right;
            
        }
        // plain vanilla and barrier in
        else {
            _dtaus.push_back(_tau_final / _Ms.front());
            _N = std::floor((_x_rs.front() - _x_ls.front()) / std::sqrt(_dtaus.front() / _alpha_temp));
            _dxs.push_back((_x_rs.front() - _x_ls.front()) / _N);
            _alphas.push_back(_dtaus.front() / (_dxs.front() * _dxs.front()));
            // get index in between value to compute
            double x = _x_ls.front(); std::size_t idx = 0;
            while (x < _target_x){
                x += _dxs.front();
                idx += 1;
            }
            _target_idx = std::make_pair(idx - 1, idx);
            
        };
    }
    // if discrete dividends
    else {
        _alphas.push_back(_alpha_temp);
        double x_bar_compute = _target_x;
        for (auto q : _q_divs){x_bar_compute += std::log(1. - q);};
        double dtau_1 = _tau_divs.front() / _Ms.front();
        _dtaus.push_back(dtau_1);
        double dx = std::sqrt(dtau_1 / _alphas.front());
        _dxs.push_back(dx);
        std::size_t N_left = std::ceil((x_bar_compute - _x_ls.front()) / dx);
        std::size_t N_right = std::ceil((_x_rs.front() - x_bar_compute) / dx);
        _target_idx = std::make_pair(N_left, N_left);
        _N = N_left + N_right;
        _x_ls[0] = x_bar_compute - N_left * _dxs.front();
        _x_rs[0] = x_bar_compute + N_right * _dxs.front();
        
        
        for (int i = 0; i < _num_divs; i++){
            _x_ls.push_back(_x_ls.back() - std::log(1. - _q_divs[i]));
            _x_rs.push_back(_x_rs.back() - std::log(1. - _q_divs[i]));
            _dtaus.push_back(_alpha_temp * _dxs.front() * _dxs.front());
            _Ms.push_back(std::ceil((_tau_divs[i+1] - _tau_divs[i]) / _dtaus.back()));
            _dtaus.back() = (_tau_divs[i+1] - _tau_divs[i]) / _Ms.back();
            _alphas.push_back(_dtaus.back() / (_dxs.back() * _dxs.back()));
            
        }
    }
};

//
void FiniteDifference::set_boundaries(){
    
    // vanilla euro call
    if (_type == vanilla && _payoff == OptionPayoff::call && _ex == OptionExercise::euro){
        _boundary_tau_0 = [=](double x)->double {return std::max(_K * std::exp(_a * x) * (std::exp(x) - 1.), 0.);};
        _boundary_x_l = [](double x, double tau)->double {return 0.;};
        _boundary_x_r = [=](double x, double tau)->double { return _K * std::exp(_a * x + _b * tau) * (std::exp(_x_rs.front() - 2. * _q * tau / (_sigma * _sigma)) - std::exp(-2. * _r * tau / (_sigma * _sigma)));};
        _early_ex_premium = [=](double x, double tau)->double { return -std::numeric_limits<double>::max(); };
    }
    // vanilla amer call
    if (_type == vanilla && _payoff == OptionPayoff::call && _ex == OptionExercise::amer){
        _boundary_tau_0 = [=](double x)->double {return std::max(_K * std::exp(_a * x) * (std::exp(x) - 1.), 0.);};
        _boundary_x_l = [](double x, double tau)->double {return 0.;};
        _boundary_x_r = [=](double x, double tau)->double { return _K * std::exp(_a * x + _b * tau) * (std::exp(x) - 1.);};
        _early_ex_premium = [=](double x, double tau)->double { return _K * std::exp(_a * x + _b * tau) * std::max(std::exp(x) - 1., 0.); };

    }
    // vanilla euro put
    if (_type == vanilla && _payoff == OptionPayoff::put && _ex == OptionExercise::euro){
        _boundary_tau_0 = [=](double x)->double {return std::max(_K * std::exp(_a * x) * (1. - std::exp(x)), 0.);};
        _boundary_x_l = [=](double x, double tau)->double { return _K * std::exp(_a * x + _b * tau) * (std::exp(-2. * _r * tau / (_sigma * _sigma)) - std::exp(_x_ls.front() - 2. * _q * tau / (_sigma * _sigma)));};
        _boundary_x_r = [](double x, double tau)->double { return 0.; };
        _early_ex_premium = [=](double x, double tau)->double { return -std::numeric_limits<double>::max(); };
        
    }
    // vanilla amer put
    if (_type == vanilla && _payoff == OptionPayoff::put && _ex == OptionExercise::amer){
        _boundary_tau_0 = [=](double x)->double {return std::max(_K * std::exp(_a * x) * (1. - std::exp(x)), 0.);};
        _boundary_x_l = [=](double x, double tau)->double { return _K * std::exp(_a * x + _b * tau) * (1. - std::exp(x));};
        _boundary_x_r = [](double x, double tau)->double { return 0.; };
        _early_ex_premium = [=](double x, double tau)->double {return _K * std::exp(_a * x + _b * tau) * std::max(1. - std::exp(x), 0.); };
    }
    // down and out euro call
    if (_type == downout && _payoff == OptionPayoff::call && _ex == OptionExercise::euro){
        _boundary_tau_0 = [=](double x)->double { return std::max(_K * std::exp(_a * x) * (std::exp(x) - 1.), 0.);};
        _boundary_x_l = [](double x, double tau)->double { return 0.;};
        _boundary_x_r = [=](double x, double tau)->double { return _K * std::exp(_a * x + _b * tau) * (std::exp(x - 2. * _q * tau / (_sigma * _sigma) ) - std::exp(-2. * _r * tau / (_sigma * _sigma)));};
        _early_ex_premium = [=](double x, double tau)->double { return -std::numeric_limits<double>::max(); };
    }
}


void FiniteDifference::build_mesh(){
    // building the mesh on x
    _x_mesh.push_back(_x_ls.front());
    for (std::size_t i = 0; i < _N; i++) {
        _x_mesh.push_back(_x_mesh.back() + _dxs.front());
    };
    // generating first layer of approximations
    _u_mesh.push_back(std::vector<double>(_N + 1));
    std::transform(_x_mesh.begin(), _x_mesh.end(), _u_mesh.back().begin(), _boundary_tau_0);
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
        new_u_mesh.push_back(std::max(
            alpha * _u_mesh.back()[pos - 1] + (1. - 2. * alpha) * _u_mesh.back()[pos] + alpha * _u_mesh.back()[pos + 1], // euro
            _early_ex_premium(_x_mesh[pos], tau)    // amer early exercise
            ));
    }
    
    // right boundary
    new_u_mesh.push_back(_boundary_x_r(xr, tau));
    
    _u_mesh.push_back(new_u_mesh);
}


void FiniteDifference::advance_impl(double alpha, double tau, double xl, double xr, const mat& L, const mat& U){
    
    std::vector<double> new_u_mesh(_N + 1);
    
    vec b(_N - 1);
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

void FiniteDifference::advance_cn_lu(double alpha, double tau, double xl, double xr, const mat& L, const mat& U, const mat& b_multiplier){
    
    std::vector<double> new_u_mesh(_N + 1);
    
    vec u(_N - 1);
    std::copy(_u_mesh.back().cbegin() + 1, _u_mesh.back().cend() - 1, u.begin());
    
    vec b(_N - 1);
    b = b_multiplier * u;
    
    // boundary conditions
    b(0) += _boundary_x_l(xl, tau) * alpha * .5 + *(_u_mesh.back().cbegin()) * alpha * .5;
    b(b.size() - 1) += _boundary_x_r(xr, tau) * alpha * .5 + *(_u_mesh.back().crbegin()) * alpha * .5;
    
    // solve linear system
    LinearSolver substituter;
    b = substituter.ForwardSub(L, b);
    b = substituter.BackwardSub(U, b);
    
    // assign to new u_mesh
    std::copy(b.cbegin(), b.cend(), new_u_mesh.begin() + 1);
    
    // Add boundary condition
    *(new_u_mesh.begin()) = _boundary_x_l(xl, tau);
    *(new_u_mesh.rbegin()) = _boundary_x_r(xr, tau);
    
    _u_mesh.push_back(new_u_mesh);
    
}


void FiniteDifference::advance_cn_sor(double alpha, double tau, double xl, double xr, const mat& A, const mat& b_multiplier){
    
    std::vector<double> new_u_mesh(_N + 1);
    
    vec u(_N - 1);
    std::copy(_u_mesh.back().cbegin() + 1, _u_mesh.back().cend() - 1, u.begin());
    
    vec b(_N - 1);
    b = b_multiplier * u;
    
    // add boundary conditions
    b(0) += _boundary_x_l(xl, tau) * alpha * .5 + *(_u_mesh.back().cbegin()) * alpha * .5;
    b(b.size() - 1) += _boundary_x_r(xr, tau) * alpha * .5 + *(_u_mesh.back().crbegin()) * alpha * .5;
    
    // yes
    vec early_ex_premium(_N - 1);
    std::copy(_x_mesh.cbegin() + 1, _x_mesh.cend() - 1, early_ex_premium.begin());
    auto find_early_ex = [this, tau](double x){ return this->_early_ex_premium(x, tau); };
    std::transform(early_ex_premium.begin(), early_ex_premium.end(), early_ex_premium.begin(), find_early_ex);
    
    
    // solve linear system
    double tolerance = std::pow(10, -6);
    double omega = 1.2;
    IterativeSolver solver(A, b, early_ex_premium);
    vec sol(b.size());
    std::tie(sol, std::ignore) = solver.SORProjected_lower(omega, StoppingCriterion::consecutive, tolerance, early_ex_premium);
    
    // assign to new u_mesh
    std::copy(sol.cbegin(), sol.cend(), new_u_mesh.begin() + 1);
    
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
        advance_expl(_alphas[sub], start_tau + _dtaus[sub] * i, _x_ls[sub], _x_rs[sub]);
    }
};

void FiniteDifference::compute_sub_domain_impl(std::size_t sub, double start_tau){
    // initial matrix
    mat A(mat::Zero(_N - 1, _N - 1));
    A(0, 0) = 1. + 2. * _alphas[sub];
    for (std::size_t i = 1; i < _N - 1; i++) {
        // Fill values
        A(i - 1, i) = -_alphas[sub];
        A(i, i - 1) = -_alphas[sub];
        A(i, i) = 1. + 2. * _alphas[sub];
    }
    
    // LU decompose A
    mat L(A.rows(), A.cols());
    mat U(A.rows(), A.cols());
    Decomposer decomposer;
    std::tie(L, U) = decomposer.lu_no_pivoting(A);
    
    
    for (int i = 1; i <= _Ms[sub]; i++){
        advance_impl(_alphas[sub], start_tau + _dtaus[sub] * i, _x_ls[sub], _x_rs[sub], L, U);
    }
};
    
void FiniteDifference::compute_sub_domain_cn_lu(std::size_t sub, double start_tau){
    mat A(mat::Zero(_N - 1, _N - 1));
    A(0, 0) = 1. + _alphas[sub];
    for (std::size_t i = 1; i < _N - 1; i++) {
        // Fill values
        A(i - 1, i) = -_alphas[sub] * .5;
        A(i, i - 1) = -_alphas[sub] * .5;
        A(i, i) = 1. + _alphas[sub];
    }
    
    // Build matrix needed for the construction of b
    mat b_multiplier(mat::Zero(A.rows(), A.cols()));   // Initialize with zero matrix
    
    b_multiplier(0, 0) = 1. - _alphas[sub];
    for (std::size_t i = 1; i < _N - 1; i++) {
        // Fill values
        b_multiplier(i - 1, i) = _alphas[sub] * .5;
        b_multiplier(i, i - 1) = _alphas[sub] * .5;
        b_multiplier(i, i) = 1. - _alphas[sub];
    }
    
    // LU decompose A
    mat L(A.rows(), A.cols());
    mat U(A.rows(), A.cols());
    Decomposer decomposer;
    std::tie(L, U) = decomposer.lu_no_pivoting(A);
    
    for (int i = 1; i <= _Ms[sub]; i++){
        advance_cn_lu(_alphas[sub], start_tau + _dtaus[sub] * i, _x_ls[sub], _x_rs[sub], L, U, b_multiplier);
    }
};
    
void FiniteDifference::compute_sub_domain_cn_sor(std::size_t sub, double start_tau){
    mat A(mat::Zero(_N - 1, _N - 1));
    A(0, 0) = 1. + _alphas[sub];
    for (std::size_t i = 1; i < _N - 1; i++) {
        // Fill values
        A(i - 1, i) = -_alphas[sub] * .5;
        A(i, i - 1) = -_alphas[sub] * .5;
        A(i, i) = 1. + _alphas[sub];
    }
    
    // Build matrix needed for the construction of b
    mat b_multiplier(mat::Zero(A.rows(), A.cols()));   // Initialize with zero matrix
    
    b_multiplier(0, 0) = 1. - _alphas[sub];
    for (std::size_t i = 1; i < _N - 1; i++) {
        // Fill values
        b_multiplier(i - 1, i) = _alphas[sub] * .5;
        b_multiplier(i, i - 1) = _alphas[sub] * .5;
        b_multiplier(i, i) = 1. - _alphas[sub];
    }
    
    for (int i = 1; i <= _Ms[sub]; i++){
        advance_cn_sor(_alphas[sub], start_tau + _dtaus[sub] * i, _x_ls[sub], _x_rs[sub], A, b_multiplier);
    }
};

void FiniteDifference::shift_domain(std::size_t sub){
    
    // shifting x and u between domains
    if (sub < _num_divs){
        // shifting _x_mesh;
        std::transform(_x_mesh.begin(), _x_mesh.end(), _x_mesh.begin(), std::bind(std::minus<double>(), std::placeholders::_1, - std::log(1. - _q_divs[sub])));
        
        // adding new u_mesh after divs
        std::vector<double> new_u_mesh(_N + 1);
        std::transform(_u_mesh.back().begin(), _u_mesh.back().end(), new_u_mesh.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, std::pow(1. - _q_divs[sub], -_a)));
        
        _u_mesh.push_back(new_u_mesh);
    }
    // shifting x back to original u to compute terminal values
    else {
        // getting total of shifts
        double shifts_sum = 0;
        for (auto q : _q_divs){ shifts_sum += std::log(1. - q);};
        std::transform(_x_mesh.begin(), _x_mesh.end(), _x_mesh.begin(), std::bind(std::plus<double>(), std::placeholders::_1,   - shifts_sum));
    }
    
}




/*
 SUPPORT
 */

double FiniteDifference::x_to_S(double x) const{
    return _K * std::exp(x);
}

double FiniteDifference::u_to_v(double x, double tau, double u) const{
    return std::exp(-_a * x -_b * tau) * u;
}

void FiniteDifference::compute_terminal_vals(){
    
    // if target is on x mesh
    if (_target_idx.first == _target_idx.second){
        
        std::size_t idx = _target_idx.first;
        
        // terminal x values
        _x_minus2 = _x_mesh[idx - 2];
        _x_minus = _x_mesh[idx - 1];
        _x_plus = _x_mesh[idx + 1];
        _x_plus2 = _x_mesh[idx + 2];
        
        // terminal u values
        _u_minus2 = _u_mesh.back()[idx - 2];
        _u_minus = _u_mesh.back()[idx - 1];
        _u_plus = _u_mesh.back()[idx + 1];
        _u_plus2 = _u_mesh.back()[idx + 2];
        
        // terminal S values
        _S_minus2 = x_to_S(_x_minus2);
        _S_minus = x_to_S(_x_minus);
        _S_plus = x_to_S(_x_plus);
        _S_plus2 = x_to_S(_x_plus2);
        
        
        // terminal V values
        _V_minus2 = u_to_v(_x_minus2, _tau_final, _u_minus2);
        _V_minus = u_to_v(_x_minus, _tau_final, _u_minus);
        _V_plus = u_to_v(_x_plus, _tau_final, _u_plus);
        _V_plus2 = u_to_v(_x_plus2, _tau_final, _u_plus2);
    }
    
    // if target is NOT on x mesh
    else {
        
        std::size_t idx_minus = _target_idx.first;
        std::size_t idx_plus = _target_idx.second;
        
        // terminal x values
        _x_minus2 = _x_mesh[idx_minus - 1];
        _x_minus = _x_mesh[idx_minus];
        _x_plus = _x_mesh[idx_plus];
        _x_plus2 = _x_mesh[idx_plus + 1];
        
        // terminal u values
        _u_minus2 = _u_mesh.back()[idx_minus - 1];
        _u_minus = _u_mesh.back()[idx_minus];
        _u_plus = _u_mesh.back()[idx_plus];
        _u_plus2 = _u_mesh.back()[idx_plus + 1];
        
        // terminal S values
        _S_minus2 = x_to_S(_x_minus2);
        _S_minus = x_to_S(_x_minus);
        _S_plus = x_to_S(_x_plus);
        _S_plus2 = x_to_S(_x_plus2);
        
        // terminal V values
        _V_minus2 = u_to_v(_x_minus2, _tau_final, _u_minus2);
        _V_minus = u_to_v(_x_minus, _tau_final, _u_minus);
        _V_plus = u_to_v(_x_plus, _tau_final, _u_plus);
        _V_plus2 = u_to_v(_x_plus2, _tau_final, _u_plus2);
    }

    
}

std::vector<double> FiniteDifference::approximate(){
    
    std::vector<double> res;
    
    std::size_t idx_minus = _target_idx.first;
    std::size_t idx_plus = _target_idx.second;
    
    // THIS NEEDS TO BE CHANGED TO REAL VALUE
    double real_val = _option.price_european()[0];
    
    // 1. linear interpolation of V
    double approx_1 = ((_S_plus - _S) * _V_minus  + (_S - _S_minus) * _V_plus) / (_S_plus - _S_minus);
    double err_1  = std::abs(approx_1 - real_val);
    
    // 2. linear interpolation of u
    double u_approx = ((_x_mesh[idx_plus] - _target_x) * _u_mesh.back()[idx_minus] + (_target_x - _x_mesh[idx_minus]) * _u_mesh.back()[idx_plus]) / _dxs.back();
    double approx_2 = std::exp(-_a * _target_x - _b * _tau_final) * u_approx;
    double err_2  = std::abs(approx_2 - real_val);
    
    // 3. RMS error
    std::vector<double> V_mesh_approx(_N + 1);
    auto calc_V_approx = [this](double x, double u) { return this->u_to_v(x, _tau_final, u); };
    std::transform(_x_mesh.begin(), _x_mesh.end(), _u_mesh.back().begin(), V_mesh_approx.begin(), calc_V_approx);
    
    std::vector<double> V_mesh_exact(_N + 1);
    // THIS NEEDS TO BE CHANGED IF ITS NOT A PUT
    auto calc_V_exact = [this](double x) {
        double S = this->x_to_S(x);
        double d1 = (std::log(S / _K) + (_r - _q + _sigma * _sigma / 2.) * _T) / (_sigma * std::sqrt(_T));
        double d2 = d1 - _sigma * std::sqrt(_T);
        return -S * std::exp(-_q * _T) * (1. - _option.phi(d1)) + _K * std::exp(-_r * _T) * (1. - _option.phi(d2));};
    std::transform(_x_mesh.cbegin(), _x_mesh.cend(), V_mesh_exact.begin(), calc_V_exact);
    
    double error_sq = 0.;
    int error_count = 0;
    auto V_mesh_approx_it = V_mesh_approx.cbegin();
    auto V_mesh_exact_it = V_mesh_exact.cbegin();
    while (V_mesh_exact_it != V_mesh_exact.cend()) {
        double V_BS = *V_mesh_exact_it;
        double V_FD = *V_mesh_approx_it;
        if (V_BS > 0.00001 * _S) {
            error_count++;
            error_sq += (V_BS - V_FD) * (V_BS - V_FD) / (V_BS * V_BS);
        }
        V_mesh_approx_it++;
        V_mesh_exact_it++;
    }
    double error_RMS = std::sqrt(error_sq / error_count);
    
    // inserts
    res.insert(res.end(), approx_1);    // approx 1 through interpolation of V
    //res.insert(res.end(), err_1);       // error approx 1
    //res.insert(res.end(), 0.);
    //res.insert(res.end(), approx_2);    // approx 2 through interpolation of u
    //res.insert(res.end(), err_2);       // error approx 2
    //res.insert(res.end(), 0.);
    //res.insert(res.end(), error_RMS);   // error RMS
    //res.insert(res.end(), 0.);
    
    
    return res;
}

std::vector<double> FiniteDifference::greeks(){

    
    std::vector<double> greeks;
    
    // if target is on a point of x_mesh
    if (_target_idx.first == _target_idx.second){
                
        double V_0 = u_to_v(_x_mesh[_target_idx.first], _tau_final, _u_mesh.back()[_target_idx.first]);
        double V_past = u_to_v(_x_mesh[_target_idx.first], _tau_final - _dtaus.back(), _u_mesh[_u_mesh.size() - 2][_target_idx.first]);
        
        double dt = 2. * _dtaus.back() / (_sigma * _sigma);
        
        double delta = (_V_plus - _V_minus) / (_S_plus - _S_minus);
        double gamma = ((_S - _S_minus) * _V_plus  - (_S_plus - _S_minus) * V_0 + (_S_plus - _S) * _V_minus) / ((_S - _S_minus) * (_S_plus - _S) * (_S_plus - _S_minus)/2.);
        double theta = (V_past - V_0) / dt;
        
        greeks.push_back(delta); greeks.push_back(gamma); greeks.push_back(theta);
        
    }
    // if target is on NOT a point of x_mesh
    else {
        
        double V_0 = ((_S_plus - _S) * _V_minus  + (_S - _S_minus) * _V_plus) / (_S_plus - _S_minus);
        double u_past_minus = _u_mesh[_u_mesh.size() - 2][_target_idx.first];
        double u_past_plus = _u_mesh[_u_mesh.size() - 2][_target_idx.second];
        double V_past_minus = u_to_v(_x_minus, _tau_final - _dtaus.back(), u_past_minus);
        double V_past_plus = u_to_v(_x_plus, _tau_final - _dtaus.back(), u_past_plus);
        double V_past = ((_S_plus - _S) * V_past_minus + (_S - _S_minus) * V_past_plus) / (_S_plus - _S_minus);
        
        double dt = 2. * _dtaus.back() / (_sigma * _sigma);
    
        double delta = (_V_plus - _V_minus) / (_S_plus - _S_minus);
        double gamma = ((_V_plus2 - _V_plus) / (_S_plus2 - _S_plus) - (_V_minus - _V_minus2) / (_S_minus - _S_minus2))
        / ((_S_plus2 + _S_plus) / 2. - (_S_minus + _S_minus2) / 2.);
        double theta = (V_past - V_0) / dt;
        
        greeks.push_back(delta); greeks.push_back(gamma); greeks.push_back(theta);
        
    }
    
    return greeks;

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



std::vector<double> FiniteDifference::price_option(const Scheme& scheme){
    
    // setting taus and xr/xl
    set_domain();
    // setting N, dx, dtau, alphas etc..
    set_discretisation();
    // setting the function at boundaries depending on option
    set_boundaries();
    // building the mesh
    build_mesh();
    
    // setting the function to compute approximations based on scheme
    std::function<void (FiniteDifference*, std::size_t, double)> scheme_func;
    switch (scheme){
        case eul_expl:
            scheme_func = &FiniteDifference::compute_sub_domain_expl;
            break;
        case eul_impl:
            scheme_func = &FiniteDifference::compute_sub_domain_impl;
            break;
        case cn_lu:
            scheme_func = &FiniteDifference::compute_sub_domain_cn_lu;
            break;
        case cn_sor:
            scheme_func = &FiniteDifference::compute_sub_domain_cn_sor;
            break;
        default:
            std::cout << "scheme not recognised" << std::endl;
    };
    
    double tau = 0;
    // computing finite differences iterively over all subdomains (divided by discrete divs)
    for (int sub = 0; sub <= _num_divs; sub++){
        // computing fd
        scheme_func(this, sub, tau);
        // shift x and u at discrete dividends (and shifts back x at tau final)
        shift_domain(sub);
        
        tau += _Ms[sub] * _dtaus[sub];
        
    }

    compute_terminal_vals();

    std::vector<double> res;
    
    // lauching approximation if target
    if (_target_idx.first != _target_idx.second){
        auto approx = approximate();
        res.insert(res.end(), approx.begin(), approx.end());
    }
    // no need to approximate as x target is already on mesh
    else {
        res.push_back(_u_mesh.back()[_target_idx.first]);
        res.push_back(u_to_v(_target_x, _tau_final, _u_mesh.back()[_target_idx.first]));
    }
    
    
    auto gs = greeks();
    res.insert(res.end(), gs.begin(), gs.end());
    
    return res;
};

void FiniteDifference::show_domain_params(){
    std::cout << "N: " << _N << std::endl;
    std::cout << "tau_final: " << _tau_final << std::endl;
    std::cout << "tau divs: "; print(_tau_divs);
    std::cout << "x_ls: "; print(_x_ls);
    std::cout  <<  "x_rs: "; print(_x_rs);
    std::cout << "dxs: "; print(_dxs);
    std::cout << "Ms: "; print(_Ms);
    std::cout << "dtaus: "; print(_dtaus);
    std::cout << "alphas: "; print(_alphas);
    std::cout << std::endl;
};

void FiniteDifference::show_grid(bool convert_to_v){
    
    std::cout << std::endl << "grid:" << std::endl;
     
    if (convert_to_v){
        
        double tau = 0;
        double x = _x_ls.front();
        
        for (int i = 0; i < _num_divs + 1; i++ ){
            for (int j = 0; j < _Ms[i]; j++){
                x = _x_ls[i];
                for (int k = 0; k < _N + 1; k++){
                    std::cout << u_to_v(x, tau, _u_mesh[j][k]) << "\t";
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

