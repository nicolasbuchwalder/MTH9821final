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

void FiniteDifference::update_params(Option opt, std::size_t M_1, double alpha_1_temp){
    std::tie(_ex, _payoff, _type, _S, _K, _T, _sigma, _r, _q, _divs, _add_params) = opt.get_params();
    
    // making sure no continuous compounding if discrete dividends
    if (_divs.empty()) _q = 0;
    
    auto num_sub_domains = _divs.size() + 1;
    _Ms.resize(num_sub_domains);
    _Ms[0] = M_1;
    _Ns.resize(num_sub_domains);
    _x_ls.resize(num_sub_domains);
    _x_rs.resize(num_sub_domains);
    _taus.resize(num_sub_domains);
    _alphas.resize(num_sub_domains);
    _alphas[0] = alpha_1_temp;
    _dxs.resize(num_sub_domains);
    _dtaus.resize(num_sub_domains);
}



void FiniteDifference::advance_expl(double alpha, double tau, double xl, double xr){
    std::vector<double> new_u_mesh;
    
    // left boundary
    new_u_mesh.push_back(_boundary_x_l_(tau, xl));
    
    // middle values
    for (std::size_t pos = 1; pos < _x_mesh.size() - 1; pos++) {
        new_u_mesh.push_back(alpha * _u_mesh.back()[pos - 1] + (1. - 2. * alpha) * _u_mesh.back()[pos] + alpha * _u_mesh.back()[pos + 1]);
    }
    
    // right boundary
    new_u_mesh.push_back(_boundary_x_r_(tau, xr));
    
    _u_mesh.push_back(new_u_mesh);
}


void FiniteDifference::advance_impl(double alpha, double tau, double xl, double xr, const mat& L, const mat& U){
    std::vector<double> new_u_mesh(_u_mesh.back().size());
    vec b(_u_mesh.back().size() - 2);
    std::copy(_u_mesh.back().cbegin() + 1, _u_mesh.back().cend() - 1, b.begin());
    
    // add boundary conditions
    b(0) += _boundary_x_l_(tau, xl) * alpha;
    b(b.size() - 1) += _boundary_x_r_(tau, xr) * alpha;
    
    // Solve linear system
    LinearSolver substituter;
    b = substituter.ForwardSub(L, b);
    b = substituter.BackwardSub(U, b);
    
    // Assign to new u_mesh
    std::copy(b.cbegin(), b.cend(), new_u_mesh.begin() + 1);
    
    // Add boundary condition
    *(new_u_mesh.begin()) = _boundary_x_l_(tau, xl);
    *(new_u_mesh.rbegin()) = _boundary_x_r_(tau, xr);
    
    _u_mesh.push_back(new_u_mesh);
}





