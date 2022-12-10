//
//  FiniteDifferencePricer.cpp
//  MTH9815_HW12_FD
//
//  Created by Nicolas Buchwalder on 28.11.22.
//

#include "FiniteDifferencePricer.hpp"
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <functional>

#include "LinearSolver.hpp"
#include "Decomposer.hpp"
#include "IterativeSolver.hpp"

//#include "LUDecomp.hpp"
//#include "LUSub.hpp"
//#include <boost/numeric/ublas/symmetric.hpp>
//#include <boost/numeric/ublas/operation.hpp>
//#include <boost/numeric/ublas/io.hpp>

#include <iostream>

double FiniteDifferencePricer::alpha_temp_ = .4;
std::size_t FiniteDifferencePricer::M_init_ = 4;

FiniteDifferencePricer::FiniteDifferencePricer(double S0, double K, double T, double sigma, double r, double q, double t_div, double B)
: S0_(S0), K_(K), T_(T), sigma_(sigma), q_(q), r_(r), t_div_(t_div), _B(B){
    
    //  1. Computational domain (part 1)
    double sigma2 = sigma * sigma;
    
    tau_div_ = (T - t_div_) * sigma2 / 2.;
    tau_final_ = T * sigma2 / 2.;
    
    x_l_ = std::log(S0 / K) + (r - sigma2 / 2.) * T - 3. * sigma * std::sqrt(T);
    x_r_ = std::log(S0 / K) + (r - sigma2 / 2.) * T + 3. * sigma * std::sqrt(T);
    
    
    // Heat equation transformation coefficients
    a_ = (r - q) / sigma2 - .5;
    b_ = (a_ + 1.) * (a_ + 1.) + 2 * q / sigma2;
    
    
}

void FiniteDifferencePricer::EuroPut_advance(double tau, double xl, double xr, double alpha, const std::vector<double>& x_mesh, std::vector<double>& u_mesh) const {
    
    std::vector<double> new_u_mesh;
    
    // Left boundary
    new_u_mesh.push_back(boundary_x_l_(tau, xl));
    
    // Middle values
    for (std::size_t pos = 1; pos < x_mesh.size() - 1; pos++) {
        new_u_mesh.push_back(alpha * u_mesh[pos - 1] + (1. - 2. * alpha) * u_mesh[pos] + alpha * u_mesh[pos + 1]);
    }
    
    // Right boundary
    new_u_mesh.push_back(boundary_x_r_(tau, xr));
    
    u_mesh = std::move(new_u_mesh);
}

void FiniteDifferencePricer::EuroPut_advance_implicit(double tau, double alpha, const std::vector<double>& x_mesh, std::vector<double>& u_mesh, const mat& L, const mat& U) const {
    
    vec b(u_mesh.size() - 2);
    std::copy(u_mesh.cbegin() + 1, u_mesh.cend() - 1, b.begin());
    
    // Add boundary conditions
    b(0) += boundary_x_l_(tau, x_l_) * alpha;
    b(b.size() - 1) += boundary_x_r_(tau, x_r_) * alpha;
    
    // Solve linear system
    LinearSolver substituter;
    b = substituter.ForwardSub(L, b);
    b = substituter.BackwardSub(U, b);
    
    // Assign to new u_mesh
    std::copy(b.cbegin(), b.cend(), u_mesh.begin() + 1);
    
    // Add boundary condition
    *(u_mesh.begin()) = boundary_x_l_(tau, x_l_);
    *(u_mesh.rbegin()) = boundary_x_r_(tau, x_r_);
    
//    this->PrintVector(u_mesh);
}

void FiniteDifferencePricer::EuroPut_advance_imex(double tau, double xl, double xr, double alpha, const std::vector<double>& x_mesh, std::vector<double>& u_mesh, const mat& L, const mat& U, const mat& b_multiplier) const {
    
    // Prepare b
    vec b(u_mesh.size() - 2);
    
    // [Some matrix] * u
    vec u(u_mesh.size() - 2);
    std::copy(u_mesh.cbegin() + 1, u_mesh.cend() - 1, u.begin());
    
//    int test_bmulsize1 = b_multiplier.size1();
//    int test_bmulsize2 = b_multiplier.size2();
//    int test_u = u.size();
//    int test_b = b.size();
    
    b = b_multiplier * u;
    
    // Boundary conditions
    b(0) += boundary_x_l_(tau, xl) * alpha * .5;
    b(b.size() - 1) += boundary_x_r_(tau, xr) * alpha * .5;
    b(0) += *(u_mesh.cbegin()) * alpha * .5;
    b(b.size() - 1) += *(u_mesh.crbegin()) * alpha * .5;
    
    // Solve linear system
    LinearSolver substituter;
    b = substituter.ForwardSub(L, b);
    b = substituter.BackwardSub(U, b);
    
    // Assign to new u_mesh
    std::copy(b.cbegin(), b.cend(), u_mesh.begin() + 1);
    
    // Add boundary condition
    *(u_mesh.begin()) = boundary_x_l_(tau, xl);
    *(u_mesh.rbegin()) = boundary_x_r_(tau, xr);
    
//    this->PrintVector(u_mesh);
}

void FiniteDifferencePricer::AmeriPut_advance_imex(double tau, double alpha, const std::vector<double>& x_mesh, std::vector<double>& u_mesh, const mat& A, const mat& b_multiplier) const {
    

    vec b(u_mesh.size() - 2);
    
    // [Some matrix] * u
    vec u(u_mesh.size() - 2);
    std::copy(u_mesh.cbegin() + 1, u_mesh.cend() - 1, u.begin());
    
//    int test_bmulsize1 = b_multiplier.size1();
//    int test_bmulsize2 = b_multiplier.size2();
//    int test_u = u.size();
//    int test_b = b.size();
    
    b = b_multiplier * u;
    
//    auto boundary = [=](double x)->double {
//        if (x < 0.) {
//            return K_ * std::exp(a_ * x + b_ * tau) * (1. - std::exp(x));
//        } else {
//            return 0.;
//        }
//    };
    
    // Boundary conditions
    b(0) += boundary_x_l_(tau, x_l_) * alpha * .5;
    b(b.size() - 1) += boundary_x_r_(tau, x_r_) * alpha * .5;
    b(0) += *(u_mesh.cbegin()) * alpha * .5;
    b(b.size() - 1) += *(u_mesh.crbegin()) * alpha * .5;
    
    double old_right_boundary = *(u_mesh.rbegin());
    
    // Add boundary condition
    *(u_mesh.begin()) = boundary_x_l_(tau, x_l_);
    *(u_mesh.rbegin()) = boundary_x_r_(tau, x_r_);
    
    auto find_early_ex = [=](double x)->double {
        if (x < 1.) {
            return 0.;//K_ * std::exp(a_ * x + b_ * tau) * (1. - std::exp(x));
        } else {
            return 0.;
        }
    };
    
    vec early_ex_premium(x_mesh.size() - 2);
    std::copy(x_mesh.cbegin() + 1, x_mesh.cend() - 1, early_ex_premium.begin());
    std::transform(early_ex_premium.begin(), early_ex_premium.end(), early_ex_premium.begin(), find_early_ex);
//    std::cout << "??" << std::endl;
//    std::cout << early_ex_premium << std::endl;
    
    
    // Solve linear system
    double tolerance = std::pow(10, -6);
    double omega = 1.2;
    // Consecutive approximation criterion
//    vec diff(vec::Ones(xold.size()));
//    vec xnew(xold.size());
//    while (diff.norm() > tolerance) {
//        xnew(0) = (1. - omega) * xold(0) + (omega * alpha) / (2. * (1. + alpha)) * (*(u_mesh.begin()) + xold(1)) + omega / (1. + alpha) * b(0);
//
//        for (std::size_t j = 1; j < xold.size() - 1; j++) {
//            xnew(j) = (1. - omega) * xold(j) + (omega * alpha) / (2. * (1. + alpha)) * (xnew(j-1) + xold(j+1)) + omega / (1. + alpha) * b(j);
//
//
//        }
//
//        std::size_t last = xold.size() - 1;
//        xnew(last) = (1. - omega) * xold(last) + (omega * alpha) / (2. * (1. + alpha)) * (xnew(last-1) + 0.) + omega / (1. + alpha) * b(last);
//
//
//        for (std::size_t k = 0; k < xnew.size(); k++) {
//            xnew(k) = std::max(xnew(k), early_ex_premium(k));
//        }
//
//        diff = xnew - xold;
//        xold = xnew;
//    }
    
    IterativeSolver solver(A, b, early_ex_premium);
    vec sol(b.size());
    std::tie(sol, std::ignore) = solver.SORProjected_lower(omega, StoppingCriterion::consecutive, tolerance, early_ex_premium);
    
    // Assign to new u_mesh
    std::copy(sol.cbegin(), sol.cend(), u_mesh.begin() + 1);

    
//    this->PrintVector(u_mesh);
}


void FiniteDifferencePricer::PrintVector(const std::vector<double>& vec) const {
    for (auto elem : vec) {
        std::cout << elem << "\t";
    }
    std::cout << std::endl;
}

std::vector<double> FiniteDifferencePricer::EuroPut(std::size_t M, const Euler& method) {
    
    // 2. Boundary conditions
    boundary_tau_0_ = [=](double x)->double {
        // K * exp(ax) * max(1 - exp(x), 0)
        if (x < 0.) {
            return K_ * std::exp(a_ * x) * (1. - std::exp(x));
        } else {
            return 0.;
        }
    };
    
    boundary_x_l_ = [=](double tau, double x)->double {
        // K * exp(ax_{left} + b\tau)(exp(-2r\tau / sigma^2) - exp(x_{left} - 2q\tau / sigma^2))
        return K_ * std::exp(a_ * x_l_ + b_ * tau) * (std::exp(-2. * r_ * tau / (sigma_ * sigma_)) - std::exp(x_l_ - 2. * q_ * tau / (sigma_ * sigma_)));
    };
    
    boundary_x_r_ = [](double tau, double x)->double {
        return 0.;
    };
    
    switch (method) {
        case Euler::Explicit:
            // Explicit (Forward Euler)
            return this->EuroPut_expl(M);
            break;
        case Euler::Implicit:
            // Implicit (Backward Euler)
            return this->EuroPut_impl(M);
            break;
        case Euler::ImEx:
            // Implicit-explicit (Crank-Nicolson)
            return this->EuroPut_imex(M);
            break;
        default:
            // Default to explicit method
            return this->EuroPut_expl(M);
            break;
    }
}

std::vector<double> FiniteDifferencePricer::EuroPut_expl(std::size_t M) const {
    
    std::vector<double> res;
    
    // 1. Computational domain (part 2)
    // M: number of time intervals on the tau-axis
    // N: number of x intervals on the x-axis
    std::size_t N;
    // Interval lengths
    double dtau, dx, alpha;
    std::tie(N, dtau, dx, alpha) = this->DomainParams(M);
    
    // 3. Finite difference scheme
    // Fill x mesh
    std::vector<double> x_mesh;
    std::size_t interval_i;
    std::tie(x_mesh, interval_i) = this->BuildMesh(N, dx);
    
    std::vector<std::vector<double>> u_meshes = this->FiniteDifference(alpha, x_mesh, M, dtau);
    std::vector<double> u_mesh = u_meshes[0];
    std::vector<double> u_mesh_prev = u_meshes[1];
    
    // 4. Pointwise convergence
    std::vector<double> approximations = this->Approximate(x_mesh, u_mesh, dx);
    res.insert(res.end(), approximations.cbegin(), approximations.cend());
    double V_approx = approximations.front();
    
    // 5. RMS error
    double error_RMS = this->RMS(x_mesh, u_mesh);
    res.push_back(error_RMS);
    
    // 6. Greeks
    std::vector<double> Greeks = this->Greeks(interval_i, x_mesh, u_mesh, u_mesh_prev, dtau, V_approx);
    res.insert(res.end(), Greeks.cbegin(), Greeks.cend());
    
    return res;
}

std::vector<double> FiniteDifferencePricer::EuroPut_impl(std::size_t M) const {
    
    std::vector<double> res;
    
    // 1. Computational domain (part 2)
    // M: number of time intervals on the tau-axis
    // N: number of x intervals on the x-axis
    std::size_t N;
    // Interval lengths
    double dtau, dx, alpha;
    std::tie(N, dtau, dx, alpha) = this->DomainParams(M);
    
    // 3. Finite difference scheme
    // Fill x mesh
    std::vector<double> x_mesh;
    std::size_t interval_i;
    std::tie(x_mesh, interval_i) = this->BuildMesh(N, dx);
    
    std::vector<std::vector<double>> u_meshes = this->FiniteDifference_impl(alpha, x_mesh, M, dtau);
    std::vector<double> u_mesh = u_meshes[0];
    std::vector<double> u_mesh_prev = u_meshes[1];
    
    // 4. Pointwise convergence
    std::vector<double> approximations = this->Approximate(x_mesh, u_mesh, dx);
    res.insert(res.end(), approximations.cbegin(), approximations.cend());
    double V_approx = approximations.front();
    
    // 5. RMS error
    double error_RMS = this->RMS(x_mesh, u_mesh);
    res.push_back(error_RMS);
    
    // 6. Greeks
    std::vector<double> Greeks = this->Greeks(interval_i, x_mesh, u_mesh, u_mesh_prev, dtau, V_approx);
    res.insert(res.end(), Greeks.cbegin(), Greeks.cend());
    
    return res;
    
}

std::vector<double> FiniteDifferencePricer::EuroPut_imex(std::size_t M) const {
    
    std::vector<double> res;
    
    // 1. Computational domain (part 2)
    // M: number of time intervals on the tau-axis
    // N: number of x intervals on the x-axis
    std::size_t N;
    // Interval lengths
    double dtau, dx, alpha;
    std::tie(N, dtau, dx, alpha) = this->DomainParams(M);
    
    // 3. Finite difference scheme
    // Fill x mesh
    std::vector<double> x_mesh;
    std::size_t interval_i;
    std::tie(x_mesh, interval_i) = this->BuildMesh(N, dx);
    
    std::vector<std::vector<double>> u_meshes = this->FiniteDifference_imex(alpha, x_mesh, M, dtau);
    std::vector<double> u_mesh = u_meshes[0];
    std::vector<double> u_mesh_prev = u_meshes[1];
    
    // 4. Pointwise convergence
    std::vector<double> approximations = this->Approximate(x_mesh, u_mesh, dx);
    res.insert(res.end(), approximations.cbegin(), approximations.cend());
    double V_approx = approximations.front();
    
    // 5. RMS error
    double error_RMS = this->RMS(x_mesh, u_mesh);
    res.push_back(error_RMS);
    
    // 6. Greeks
    std::vector<double> Greeks = this->Greeks(interval_i, x_mesh, u_mesh, u_mesh_prev, dtau, V_approx);
    res.insert(res.end(), Greeks.cbegin(), Greeks.cend());
    
    return res;
    
}



std::tuple<std::size_t, double, double, double> FiniteDifferencePricer::DomainParams(std::size_t M) const {

    // More finite difference parameters
    double dtau = tau_final_ / M;
    // N: number of x intervals on the x-axis
    std::size_t N = std::floor((x_r_ - x_l_) / std::sqrt(dtau / alpha_temp_));
    //    std::cout << N << std::endl;
    double dx = (x_r_ - x_l_) / N;
    double alpha = dtau / (dx * dx);
    
    return std::make_tuple(N, dtau, dx, alpha);
}
std::tuple<size_t, double, double, double, double> FiniteDifferencePricer::DomainParams_barrier(std::size_t M) {
    double x_compute = std::log(S0_/K_);
    x_l_ = std::log(_B/K_);
    double d_tau = tau_final_ / M;
    double dx = std::sqrt(d_tau / alpha_temp_);
    std::size_t N_left = std::floor((x_compute - x_l_) / dx);
    dx = (x_compute - x_l_) / N_left;
    double alpha = d_tau / (dx * dx);
    std::size_t N_right = std::ceil((x_r_ - x_compute) / dx);
    std::size_t N = N_left + N_right;
    
    //std::cout << "M:  \talpha    \tx_left\t\tx_right\tN\tdelta_x\tdelta_tau" << std::endl;
    //std::cout << M << "\t\t" << alpha << "\t" << x_l_ << "\t" << x_r_ << "\t" << N << "\t"<< dx << "\t" << d_tau << std::endl;
    
    return std::make_tuple(N, dx, d_tau, alpha, x_compute);
    
}

std::tuple<std::size_t, double, double, double, double, double, double, double> FiniteDifferencePricer::DomainParams_discrete_divs(std::size_t M_1) {

    double x_compute = std::log(S0_/K_);
    double x_bar_compute = x_compute + std::log(1 - q_);
    double dtau_1 = tau_div_ / M_1;
    double dx = std::sqrt(dtau_1 / alpha_temp_);
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
    
    //std::cout << "M_1:\tM_2\t  alpha_2\tN\tx_left\t\tx_right\t\tx_left,new\tx_right,new\t\ttau_div\t\tdelta_tau_1\tdelta_tau_2\tdelta_x" << std::endl;
    //std::cout << M_1 << "\t\t" << M_2 << "\t  " << alpha_2 << "\t" << N << "\t" << x_l_ << "\t" << x_r_ << "\t" << x_l_new_ << "\t\t" << x_r_new_ << "\t" << tau_div_ << "\t" << dtau_1 << "\t" << dtau_2 << "\t" << dx << std::endl;
    
    idx_compute = N_left;
    return std::make_tuple(N, M_2, dx, dtau_1, dtau_2, alpha_2, x_compute, x_bar_compute);
};


std::pair<std::vector<double>, std::size_t> FiniteDifferencePricer::BuildMesh(std::size_t N, double dx) const {
    
    // Fill x mesh
    std::vector<double> x_mesh({x_l_});
    for (std::size_t i = 0; i < N; i++) {
        x_mesh.push_back(x_mesh.back() + dx);
    }
    
    // Find actual x on x_mesh
    double x_compute = std::log(S0_ / K_);
    auto x_large_it = std::upper_bound(x_mesh.cbegin(), x_mesh.cend(), x_compute);
    std::size_t interval_i = std::distance(x_mesh.cbegin(), x_large_it) - 1;
    
    return std::make_pair(x_mesh, interval_i);
};

std::vector<double> FiniteDifferencePricer::BuildMesh_discrete_divs(std::size_t N, double dx) const{
    std::vector<double> x_mesh({x_l_});
    for (std::size_t i = 0; i < N; i++) {
        x_mesh.push_back(x_mesh.back() + dx);
        
    };
    return x_mesh;
}


std::vector<std::vector<double>> FiniteDifferencePricer::FiniteDifference(double alpha, const std::vector<double>& x_mesh, std::size_t M, double dtau) const {
    
    // Fill u mesh with boundary condition
    std::vector<double> u_mesh(x_mesh.size());
    std::transform(x_mesh.cbegin(), x_mesh.cend(), u_mesh.begin(), boundary_tau_0_);
    
    //PrintVector(u_mesh);
    // Advance M-1 times
    for (std::size_t i = 1; i < M; i++) {
        double curr_tau = dtau * i;
        this->EuroPut_advance(curr_tau, x_l_, x_r_, alpha, x_mesh, u_mesh);
        //PrintVector(u_mesh);
    }
    
    // Record the (M-1)th u_mesh for theta calculation
    std::vector<double> u_mesh_prev = u_mesh;
    this->EuroPut_advance(tau_final_, x_l_, x_r_, alpha, x_mesh, u_mesh);
    //PrintVector(u_mesh);
    
    return std::vector<std::vector<double>>({u_mesh, u_mesh_prev});
}

std::vector<std::vector<double>> FiniteDifferencePricer::FiniteDifference_impl(double alpha, const std::vector<double>& x_mesh, std::size_t M, double dtau) const {
    
    // Fill u mesh with boundary condition
    std::vector<double> u_mesh(x_mesh.size());
    std::transform(x_mesh.cbegin(), x_mesh.cend(), u_mesh.begin(), boundary_tau_0_);
    *(u_mesh.begin()) = boundary_x_l_(0., x_l_);
    *(u_mesh.rbegin()) = boundary_x_r_(0., x_r_);
//    this->PrintVector(u_mesh);
    //PrintVector(u_mesh);
    // Build matrix A
    mat A(mat::Zero(x_mesh.size() - 2, x_mesh.size() - 2));   // Initialize with zero matrix
    A(0, 0) = 1. + 2. * alpha;
    for (std::size_t i = 1; i < x_mesh.size() - 2; i++) {
        // Fill values
        A(i - 1, i) = -alpha;
        A(i, i - 1) = -alpha;
        A(i, i) = 1. + 2. * alpha;
    }
    
    // LU decompose A
    mat L(A.rows(), A.cols());
    mat U(A.rows(), A.cols());
    
    Decomposer decomposer;
    std::tie(L, U) = decomposer.lu_no_pivoting(A);
    
    // Advance M-1 times
    for (std::size_t i = 1; i < M; i++) {
        double curr_tau = dtau * i;
        this->EuroPut_advance_implicit(curr_tau, alpha, x_mesh, u_mesh, L, U);
        //PrintVector(u_mesh);
    }
    
    // Record the (M-1)th u_mesh for theta calculation
    std::vector<double> u_mesh_prev = u_mesh;
    this->EuroPut_advance_implicit(tau_final_, alpha, x_mesh, u_mesh, L, U);
    //PrintVector(u_mesh);
    
    return std::vector<std::vector<double>>({u_mesh, u_mesh_prev});
}

std::vector<std::vector<double>> FiniteDifferencePricer::FiniteDifference_imex(double alpha, const std::vector<double>& x_mesh, std::size_t M, double dtau) const {
    
    // Fill u mesh with boundary condition
    std::vector<double> u_mesh(x_mesh.size());
    std::transform(x_mesh.cbegin(), x_mesh.cend(), u_mesh.begin(), boundary_tau_0_);
//    this->PrintVector(u_mesh);
    
    mat A(mat::Zero(x_mesh.size() - 2, x_mesh.size() - 2));   // Initialize with zero matrix
    A(0, 0) = 1. + alpha;
    for (std::size_t i = 1; i < x_mesh.size() - 2; i++) {
        // Fill values
        A(i - 1, i) = -alpha * .5;
        A(i, i - 1) = -alpha * .5;
        A(i, i) = 1. + alpha;
    }
    
    // Build matrix needed for the construction of b
    mat b_multiplier(mat::Zero(A.rows(), A.cols()));   // Initialize with zero matrix
    
    b_multiplier(0, 0) = 1. - alpha;
    for (std::size_t i = 1; i < x_mesh.size() - 2; i++) {
        // Fill values
        b_multiplier(i - 1, i) = alpha * .5;
        b_multiplier(i, i - 1) = alpha * .5;
        b_multiplier(i, i) = 1. - alpha;
    }
    
    // LU decompose A
    mat L(A.rows(), A.cols());
    mat U(A.rows(), A.cols());
    
    Decomposer decomposer;
    std::tie(L, U) = decomposer.lu_no_pivoting(A);
    
    
    // Advance M-1 times
    for (std::size_t i = 1; i < M; i++) {
        double curr_tau = dtau * i;
        this->EuroPut_advance_imex(curr_tau, 0, 0, alpha, x_mesh, u_mesh, L, U, b_multiplier);
    }
    
    // Record the (M-1)th u_mesh for theta calculation
    std::vector<double> u_mesh_prev = u_mesh;
    this->EuroPut_advance_imex(tau_final_, 0, 0, alpha, x_mesh, u_mesh, L, U, b_multiplier);
    
    return std::vector<std::vector<double>>({u_mesh, u_mesh_prev});
}

std::vector<std::vector<double>> FiniteDifferencePricer::FiniteDifference_imex_AmeriPut(double alpha, const std::vector<double>& x_mesh, std::size_t M, double dtau) const {
    
    // Fill u mesh with boundary condition
    std::vector<double> u_mesh(x_mesh.size());
    std::transform(x_mesh.cbegin(), x_mesh.cend(), u_mesh.begin(), boundary_tau_0_);
    
//    this->PrintVector(u_mesh);
    
    mat A(mat::Zero(x_mesh.size() - 2, x_mesh.size() - 2));   // Initialize with zero matrix
    A(0, 0) = 1. + alpha;
    for (std::size_t i = 1; i < x_mesh.size() - 2; i++) {
        // Fill values
        A(i - 1, i) = -alpha * .5;
        A(i, i - 1) = -alpha * .5;
        A(i, i) = 1. + alpha;
    }
    
    // Build matrix needed for the construction of b
    mat b_multiplier(mat::Zero(A.rows(), A.cols()));   // Initialize with zero matrix
    
    b_multiplier(0, 0) = 1. - alpha;
    for (std::size_t i = 1; i < x_mesh.size() - 2; i++) {
        // Fill values
        b_multiplier(i - 1, i) = alpha * .5;
        b_multiplier(i, i - 1) = alpha * .5;
        b_multiplier(i, i) = 1. - alpha;
    }
    
//    // LU decompose A
//    mat L(A.rows(), A.cols());
//    mat U(A.rows(), A.cols());
//
//    Decomposer decomposer;
//    std::tie(L, U) = decomposer.lu_no_pivoting(A);
    
    
    // Advance M-1 times
    for (std::size_t i = 1; i < M; i++) {
        double curr_tau = dtau * i;
        this->AmeriPut_advance_imex(curr_tau, alpha, x_mesh, u_mesh, A, b_multiplier);
    }
    
    // Record the (M-1)th u_mesh for theta calculation
    std::vector<double> u_mesh_prev = u_mesh;
    this->AmeriPut_advance_imex(tau_final_, alpha, x_mesh, u_mesh, A, b_multiplier);
    
    return std::vector<std::vector<double>>({u_mesh, u_mesh_prev});
}

std::vector<std::vector<double>> FiniteDifferencePricer::FiniteDifference_expl_discr(const std::vector<double>& x_mesh, std::size_t M_1, double dtau_1, std::size_t M_2, double dtau_2, double alpha_2) const {

    // Fill u mesh with boundary condition
    std::vector<double> u_mesh(x_mesh.size());
    std::transform(x_mesh.cbegin(), x_mesh.cend(), u_mesh.begin(), boundary_tau_0_);
    //PrintVector(u_mesh);
//    std::vector<double> temp_mesh(u_mesh.size());
//    auto back_to_v = [=](double x, double u){ return std::exp(-a_ * x) * u;};
//    std::transform(x_mesh.begin(), x_mesh.end(), u_mesh.begin(), temp_mesh.begin(), back_to_v);
//    PrintVector(temp_mesh);
    
    
    // Advance M_1 times
    for (std::size_t i = 1; i <= M_1; i++) {
        double curr_tau = dtau_1 * i;
        this->EuroPut_advance(curr_tau, x_l_, x_r_, alpha_temp_, x_mesh, u_mesh);
        //PrintVector(u_mesh);
    }
    
    
    std::transform(u_mesh.begin(), u_mesh.end(), u_mesh.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, std::pow(1-q_, -a_)));
    //PrintVector(u_mesh);
    
    // Advance M_2-1 times
    for (std::size_t i = 1; i < M_2; i++) {
        double curr_tau = tau_div_ + dtau_2 * i;
        this->EuroPut_advance(curr_tau, x_l_new_, x_r_new_, alpha_2, x_mesh, u_mesh);
        //PrintVector(u_mesh);
    }
    // Record the (M-1)th u_mesh for theta calculation
    std::vector<double> u_mesh_prev = u_mesh;
    this->EuroPut_advance(tau_final_, x_l_new_, x_r_new_, alpha_2, x_mesh, u_mesh);
    
    //PrintVector(u_mesh);

    return std::vector<std::vector<double>>({u_mesh, u_mesh_prev});
}

std::vector<std::vector<double>> FiniteDifferencePricer::FiniteDifference_discrete_imex(const std::vector<double>& x_mesh, std::size_t M_1, double dtau_1, std::size_t M_2, double dtau_2, double alpha_2) const {
    
    // Fill u mesh with boundary condition
    std::vector<double> u_mesh(x_mesh.size());
    std::transform(x_mesh.cbegin(), x_mesh.cend(), u_mesh.begin(), boundary_tau_0_);
    
//    this->PrintVector(u_mesh);
    
    mat A(mat::Zero(x_mesh.size() - 2, x_mesh.size() - 2));   // Initialize with zero matrix
    A(0, 0) = 1. + alpha_temp_;
    for (std::size_t i = 1; i < x_mesh.size() - 2; i++) {
        // Fill values
        A(i - 1, i) = -alpha_temp_ * .5;
        A(i, i - 1) = -alpha_temp_ * .5;
        A(i, i) = 1. + alpha_temp_;
    }
    
    // Build matrix needed for the construction of b
    mat b_multiplier(mat::Zero(A.rows(), A.cols()));   // Initialize with zero matrix
    
    b_multiplier(0, 0) = 1. - alpha_temp_;
    for (std::size_t i = 1; i < x_mesh.size() - 2; i++) {
        // Fill values
        b_multiplier(i - 1, i) = alpha_temp_ * .5;
        b_multiplier(i, i - 1) = alpha_temp_ * .5;
        b_multiplier(i, i) = 1. - alpha_temp_;
    }
    
    // LU decompose A
    mat L(A.rows(), A.cols());
    mat U(A.rows(), A.cols());
    
    Decomposer decomposer;
    std::tie(L, U) = decomposer.lu_no_pivoting(A);
    
    \
    // Advance M_1 times
    for (std::size_t i = 1; i <= M_1; i++) {
        double curr_tau = dtau_1 * i;
        this->EuroPut_advance_imex(curr_tau, x_l_, x_r_, alpha_temp_, x_mesh, u_mesh, L, U, b_multiplier);
    }
    
    std::transform(u_mesh.begin(), u_mesh.end(), u_mesh.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, std::pow(1-q_, -a_)));
    
    
    A(0, 0) = 1. + alpha_2;
    for (std::size_t i = 1; i < x_mesh.size() - 2; i++) {
        // Fill values
        A(i - 1, i) = -alpha_2 * .5;
        A(i, i - 1) = -alpha_2 * .5;
        A(i, i) = 1. + alpha_2;
    }
    
    // Build matrix needed for the construction of b
    
    b_multiplier(0, 0) = 1. - alpha_2;
    for (std::size_t i = 1; i < x_mesh.size() - 2; i++) {
        // Fill values
        b_multiplier(i - 1, i) = alpha_2 * .5;
        b_multiplier(i, i - 1) = alpha_2 * .5;
        b_multiplier(i, i) = 1. - alpha_2;
    }
    
    std::tie(L, U) = decomposer.lu_no_pivoting(A);
    
    
    // Advance M_2-1 times
    for (std::size_t i = 1; i < M_2; i++) {
        double curr_tau = tau_div_ + dtau_2 * i;
        this->EuroPut_advance_imex(curr_tau, x_l_new_, x_r_new_, alpha_2, x_mesh, u_mesh, L, U, b_multiplier);
    }
    
    // Record the (M-1)th u_mesh for theta calculation
    std::vector<double> u_mesh_prev = u_mesh;
    this->EuroPut_advance_imex(tau_final_, 0, 0, alpha_2, x_mesh, u_mesh, L, U, b_multiplier);
    
    return std::vector<std::vector<double>>({u_mesh, u_mesh_prev});
}

std::vector<double> FiniteDifferencePricer::Approximate(const std::vector<double>& x_mesh, const std::vector<double>& u_mesh, double dx) const {
    
    std::vector<double> approximations;
    
    // Find the interval containing x_compute
    double x_compute = std::log(S0_ / K_);
    auto x_large_it = std::upper_bound(x_mesh.cbegin(), x_mesh.cend(), x_compute);
    std::size_t interval_i = std::distance(x_mesh.cbegin(), x_large_it) - 1;
    
    // Approximation method 1
    double S_small = K_ * std::exp(x_mesh[interval_i]);
    double S_large = K_ * std::exp(x_mesh[interval_i + 1]);
    double V_small = std::exp(-a_ * (x_mesh[interval_i]) - b_ * tau_final_) * u_mesh[interval_i];
    double V_large = std::exp(-a_ * (x_mesh[interval_i + 1]) - b_ * tau_final_) * u_mesh[interval_i + 1];
    double V_approx = ((S_large - S0_) * V_small + (S0_ - S_small) * V_large) / (S_large - S_small);
    approximations.push_back(V_approx);
    
    // Approximation method 2
    double u_approx = ((x_mesh[interval_i + 1] - x_compute) * u_mesh[interval_i] + (x_compute - x_mesh[interval_i]) * u_mesh[interval_i + 1]) / (dx);
    double V_approx_2 = std::exp(-a_ * x_compute - b_ * tau_final_) * u_approx;
    approximations.push_back(V_approx_2);
    
    return approximations;
}

double FiniteDifferencePricer::RMS(const std::vector<double>& x_mesh, const std::vector<double>& u_mesh) const {
    
    // Calculate vector of FD and BS option values
    auto V_mesh_approx_gen = [&](double x, double u)->double {
        return std::exp(-a_ * x - b_ * tau_final_) * u;
    };
    auto V_mesh_exact_gen = [&](double x)->double {
        return 4.219638;
    };
    
    std::vector<double> V_mesh_approx(x_mesh.size());
    std::vector<double> V_mesh_exact(x_mesh.size());
    std::transform(x_mesh.cbegin(), x_mesh.cend(), u_mesh.cbegin(), V_mesh_approx.begin(), V_mesh_approx_gen);
    std::transform(x_mesh.cbegin(), x_mesh.cend(), V_mesh_exact.begin(), V_mesh_exact_gen);
    
    // Find RMS
    double error_sq = 0.;
    int error_count = 0;
    auto V_mesh_approx_it = V_mesh_approx.cbegin();
    auto V_mesh_exact_it = V_mesh_exact.cbegin();
    while (V_mesh_exact_it != V_mesh_exact.cend()) {
        double V_BS = *V_mesh_exact_it;
        double V_FD = *V_mesh_approx_it;
        if (V_BS > 0.00001 * S0_) {
            error_count++;
            error_sq += (V_BS - V_FD) * (V_BS - V_FD) / (V_BS * V_BS);
        }
        V_mesh_approx_it++;
        V_mesh_exact_it++;
    }
    double error_RMS = std::sqrt(error_sq / error_count);
    
    return error_RMS;
}

std::vector<double> FiniteDifferencePricer::Greeks(std::size_t interval_i, const std::vector<double>& x_mesh, const std::vector<double>& u_mesh, std::vector<double>& u_mesh_prev, double dtau, double V_approx) const {
    
    // DELTA & GAMMA
    // Get S and V of interest
    double S_smaller = K_ * std::exp(x_mesh[interval_i - 1]);
    double S_small = K_ * std::exp(x_mesh[interval_i]);
    double S_large = K_ * std::exp(x_mesh[interval_i + 1]);
    double S_larger = K_ * std::exp(x_mesh[interval_i + 2]);
    
    double V_smaller = std::exp(-a_ * (x_mesh[interval_i - 1]) - b_ * tau_final_) * u_mesh[interval_i - 1];
    double V_small = std::exp(-a_ * (x_mesh[interval_i]) - b_ * tau_final_) * u_mesh[interval_i];
    double V_large = std::exp(-a_ * (x_mesh[interval_i + 1]) - b_ * tau_final_) * u_mesh[interval_i + 1];
    double V_larger = std::exp(-a_ * (x_mesh[interval_i + 2]) - b_ * tau_final_) * u_mesh[interval_i + 2];
    
    // Find delta and gamma
    double delta = (V_large - V_smaller) / (S_large - S_smaller);
    //double gamma = ((V_larger - V_large) / (S_larger - S_large) - (V_small - V_smaller) / (S_small - S_smaller)) / (((S_larger + S_large) / 2.) - ((S_small + S_smaller) / 2.));
    double gamma = ((S_small - S_smaller) * V_large  - (S_large - S_smaller) * V_small + (S_large - S_small) * V_smaller) / ((S_small - S_smaller) * (S_large - S_small) * (S_large - S_smaller)/2.);
    // THETA
    // Get dt from dtau
    double dt = 2. * dtau / (sigma_ * sigma_);
    
    // Get V at t = dt
    double V_small_prev = std::exp(-a_ * (x_mesh[interval_i]) - b_ * (tau_final_ - dtau)) * u_mesh_prev[interval_i];
    double V_large_prev = std::exp(-a_ * (x_mesh[interval_i + 1]) - b_ * (tau_final_ - dtau)) * u_mesh_prev[interval_i + 1];
    double V_approx_prev = ((S_large - S0_) * V_small_prev + (S0_ - S_small) * V_large_prev) / (S_large - S_small);
    double theta = (V_approx_prev - V_approx) / dt;
    
    return std::vector<double>({delta, gamma, theta});
}


std::vector<double> FiniteDifferencePricer::AmericanPut(std::size_t M, const Euler& method) {
    
    // 2. Boundary conditions
    boundary_tau_0_ = [=](double x)->double {
        // K * exp(ax) * max(1 - exp(x), 0)
        if (x < 0.) {
            return K_ * std::exp(a_ * x) * (1. - std::exp(x));
        } else {
            return 0.;
        }
    };
    
    boundary_x_l_ = [=](double tau, double x)->double {
        // K * exp(ax_{left} + b\tau)(exp(-2r\tau / sigma^2) - exp(x_{left} - 2q\tau / sigma^2))
        return K_ * std::exp(a_ * x_l_ + b_ * tau) * (1. - std::exp(x_l_));
    };
    
    boundary_x_r_ = [](double tau, double x)->double {
        return 0.;
    };
    
    switch (method) {
        case Euler::Explicit:
            return this->AmericanPut_expl(M);
            break;
        case Euler::ImEx:
            return this->AmericanPut_imex(M);
            break;
        default:
            return this->AmericanPut_expl(M);
            break;
    }
}

std::vector<double> FiniteDifferencePricer::AmericanPut_expl(std::size_t M) const {
    
    std::vector<double> res;
    
    // 1. Computational domain (part 2)
    // M: number of time intervals on the tau-axis
    // N: number of x intervals on the x-axis
    std::size_t N;
    // Interval lengths
    double dtau, dx, alpha;
    std::tie(N, dtau, dx, alpha) = this->DomainParams(M);
    
    // 3. Finite difference scheme
    // Fill x mesh
    std::vector<double> x_mesh;
    std::size_t interval_i;
    std::tie(x_mesh, interval_i) = this->BuildMesh(N, dx);
    
    std::vector<std::vector<double>> u_meshes = this->FiniteDifference_AmeriPut(alpha, x_mesh, M, dtau);
    std::vector<double> u_mesh = u_meshes[0];
    std::vector<double> u_mesh_prev = u_meshes[1];
    
    // 4. Pointwise convergence
    std::vector<double> approximations = this->Approximate(x_mesh, u_mesh, dx);
    res.insert(res.end(), approximations.cbegin(), approximations.cend());
    double V_approx = approximations.front();
    
    // 6. Greeks
    std::vector<double> Greeks = this->Greeks(interval_i, x_mesh, u_mesh, u_mesh_prev, dtau, V_approx);
    res.insert(res.end(), Greeks.cbegin(), Greeks.cend());
    
    // 7. Variance reduction
    res.push_back(this->VarianceReduction_AmeriPut(V_approx, M));
    
    return res;
    
}

std::vector<double> FiniteDifferencePricer::AmericanPut_imex(std::size_t M) const {
    
    std::vector<double> res;
    
    // 1. Computational domain (part 2)
    // M: number of time intervals on the tau-axis
    // N: number of x intervals on the x-axis
    std::size_t N;
    // Interval lengths
    double dtau, dx, alpha;
    std::tie(N, dtau, dx, alpha) = this->DomainParams(M);
    
    // 3. Finite difference scheme
    // Fill x mesh
    std::vector<double> x_mesh;
    std::size_t interval_i;
    std::tie(x_mesh, interval_i) = this->BuildMesh(N, dx);
    
    std::vector<std::vector<double>> u_meshes = this->FiniteDifference_imex_AmeriPut(alpha, x_mesh, M, dtau);
    std::vector<double> u_mesh = u_meshes[0];
    std::vector<double> u_mesh_prev = u_meshes[1];
    
    // 4. Pointwise convergence
    std::vector<double> approximations = this->Approximate(x_mesh, u_mesh, dx);
    res.insert(res.end(), approximations.cbegin(), approximations.cend());
    double V_approx = approximations.front();
    
    // 6. Greeks
    std::vector<double> Greeks = this->Greeks(interval_i, x_mesh, u_mesh, u_mesh_prev, dtau, V_approx);
    res.insert(res.end(), Greeks.cbegin(), Greeks.cend());
    
    // 7. Variance reduction
    res.push_back(this->VarianceReduction_imex_AmeriPut(V_approx, M));
    
    
    return res;
    
}

std::vector<std::vector<double>> FiniteDifferencePricer::FiniteDifference_AmeriPut(double alpha, const std::vector<double>& x_mesh, std::size_t M, double dtau) const {
    
    // Fill u mesh with boundary condition
    std::vector<double> u_mesh(x_mesh.size());
    std::transform(x_mesh.cbegin(), x_mesh.cend(), u_mesh.begin(), boundary_tau_0_);
    
    // Advance M-1 times
    for (std::size_t i = 1; i < M; i++) {
        double curr_tau = dtau * i;
        this->AmeriPut_advance(curr_tau, alpha, x_mesh, u_mesh);
    }
    
    // Record the (M-1)th u_mesh for theta calculation
    std::vector<double> u_mesh_prev = u_mesh;
    this->AmeriPut_advance(tau_final_, alpha, x_mesh, u_mesh);
    
    return std::vector<std::vector<double>>({u_mesh, u_mesh_prev});
}

void FiniteDifferencePricer::AmeriPut_advance(double tau, double alpha, const std::vector<double>& x_mesh, std::vector<double>& u_mesh) const {
    
    std::vector<double> new_u_mesh;
    
    // Left boundary
    new_u_mesh.push_back(boundary_x_l_(tau, 0.));
    
    // Middle values
    for (std::size_t pos = 1; pos < x_mesh.size() - 1; pos++) {
        // Get the corresponding European option's value
        double euro_val = alpha * u_mesh[pos - 1] + (1. - 2. * alpha) * u_mesh[pos] + alpha * u_mesh[pos + 1];
        
        // Find early exercise
        double early_ex_premium = 0.;
        
        if (x_mesh[pos] < 0.) {
            early_ex_premium = K_ * std::exp(a_ * x_mesh[pos] + b_ * tau) * (1. - std::exp(x_mesh[pos]));
        }
        
        // Compare and add to mesh
        new_u_mesh.push_back(std::max(euro_val, early_ex_premium));
    }
    
    // Right boundary
    new_u_mesh.push_back(boundary_x_r_(tau, 0.));
    
    u_mesh = std::move(new_u_mesh);
}


double FiniteDifferencePricer::VarianceReduction_AmeriPut(double V_approx, std::size_t M) const {
    
    // Get corresponding European option value for BS and FD
    FiniteDifferencePricer FDPricer(S0_, K_, T_, sigma_, r_, q_, t_div_, 0);
    //EuropeanOption BSPricer(0., S0_, K_, T_, sigma_, r_, q_);
    
    //double FDEuro = FDPricer.EuroPut(M).front();
    //double BSEuro = BSPricer.Put();
    
    // Adjust approximation by the pointwise difference
    //return V_approx + BSEuro - FDEuro;
    return 0;
    
}

double FiniteDifferencePricer::VarianceReduction_imex_AmeriPut(double V_approx, std::size_t M) const {
    
    // Get corresponding European option value for BS and FD
    //FiniteDifferencePricer FDPricer(S0_, K_, T_, sigma_, r_, q_, t_div_, 0);
    //EuropeanOption BSPricer(0., S0_, K_, T_, sigma_, r_, q_);
    
    //double FDEuro = FDPricer.EuroPut(M, Euler::ImEx).front();
    //double BSEuro = BSPricer.Put();
    
    // Adjust approximation by the pointwise difference
    //return V_approx + BSEuro - FDEuro;
    return 0;
    
}


std::array<std::vector<double>, 2> FiniteDifferencePricer::AmericanPut_EarlyExDomain(std::size_t M) {
    
    // 2. Boundary conditions
    boundary_tau_0_ = [=](double x)->double {
        // K * exp(ax) * max(1 - exp(x), 0)
        if (x < 0.) {
            return K_ * std::exp(a_ * x) * (1. - std::exp(x));
        } else {
            return 0.;
        }
    };
    
    boundary_x_l_ = [=](double tau, double x)->double {
        // K * exp(ax_{left} + b\tau)(exp(-2r\tau / sigma^2) - exp(x_{left} - 2q\tau / sigma^2))
        return K_ * std::exp(a_ * x_l_ + b_ * tau) * (std::exp(-2. * r_ * tau / (sigma_ * sigma_)) - std::exp(x_l_ - 2. * q_ * tau / (sigma_ * sigma_)));
    };
    
    boundary_x_r_ = [](double tau, double x)->double {
        return 0.;
    };
    
    // 1. Computational domain (part 2)
    // M: number of time intervals on the tau-axis
    // N: number of x intervals on the x-axis
    std::size_t N;
    // Interval lengths
    double dtau, dx, alpha;
    std::tie(N, dtau, dx, alpha) = this->DomainParams(M);
    
    // 3. Finite difference scheme
    // Fill x mesh
    std::vector<double> x_mesh;
    std::size_t interval_i;
    std::tie(x_mesh, interval_i) = this->BuildMesh(N, dx);
    
    // MODIFIED FINITE DIFFERENCE
    // Fill u mesh with boundary condition
    std::vector<double> u_mesh(x_mesh.size());
    std::transform(x_mesh.cbegin(), x_mesh.cend(), u_mesh.begin(), boundary_tau_0_);
//    std::cout << N << std::endl << std::endl;
    // Advance M times while recording early exercise position
    std::vector<double> Sopt;
    for (std::size_t i = 1; i <= M; i++) {
        // Record maximum exercise position
        std::size_t Nopt = 0;
        
        double tau = dtau * i;
        
        std::vector<double> new_u_mesh;
        
        // Left boundary
        new_u_mesh.push_back(boundary_x_l_(tau, 0.));
        
        // Middle values
        for (std::size_t pos = 1; pos < x_mesh.size() - 1; pos++) {
            // Get the corresponding European option's value
            double euro_val = alpha * u_mesh[pos - 1] + (1. - 2. * alpha) * u_mesh[pos] + alpha * u_mesh[pos + 1];
            
            // Find early exercise
            double early_ex_premium = 0.;
            
            if (x_mesh[pos] < 0.) {
                early_ex_premium = K_ * std::exp(a_ * x_mesh[pos] + b_ * tau) * (1. - std::exp(x_mesh[pos]));
            }
            
            // Compare and add to mesh
            if (early_ex_premium >= euro_val) {
                new_u_mesh.push_back(early_ex_premium);
                if (early_ex_premium > 0) Nopt = pos;
            } else {
                new_u_mesh.push_back(euro_val);
            }
        }
        double S_small = K_ * std::exp(x_mesh[Nopt]);
        double S_large = K_ * std::exp(x_mesh[Nopt + 1]);
        Sopt.push_back((S_small + S_large) / 2.);
        
        // Right boundary
        new_u_mesh.push_back(boundary_x_r_(tau, 0.));
        
        u_mesh = std::move(new_u_mesh);
    }
    
    // Get array of t
    std::vector<double> t_mesh(16);
    std::iota(t_mesh.begin(), t_mesh.end(), 1.);
    auto convert_to_t = [&](double m)->double {
        return T_ - (2. * m * dtau) / (sigma_ * sigma_);
    };
    std::transform(t_mesh.begin(), t_mesh.end(), t_mesh.begin(), convert_to_t);
    
//    this->PrintVector(x_mesh);
    return std::array<std::vector<double>, 2>({t_mesh, Sopt});
}

std::vector<double> FiniteDifferencePricer::EuroCallDiscreteQ(std::size_t M_1, const Euler& method){
    
    // 2. Boundary conditions
    boundary_tau_0_ = [=](double x)->double {
        // K * exp(ax) * max(1 - exp(x), 0)
        if (x > 0.) {
            return K_ * std::exp(a_ * x) * (std::exp(x) - 1.);
        } else {
            return 0.;
        }
    };
    
    boundary_x_l_ = [](double tau, double x)->double {
        return 0.;
    };
    
    boundary_x_r_ = [=](double tau, double x)->double {
        if (x == x_r_){
            return K_ * std::exp(a_ * x_r_ + b_ * tau) * (std::exp(x_r_) - std::exp(-2. * r_ * tau / (sigma_ * sigma_)));
        }
        // K * exp(ax_{left} + b\tau)(exp(-2r\tau / sigma^2) - exp(x_{left} - 2q\tau / sigma^2))
        else {
            return K_ * std::exp(a_ * x_r_new_ + b_ * tau) * (std::exp(x_r_) - std::exp(-2. * r_ * tau / (sigma_ * sigma_)));
        }
    };
    
    switch (method) {
        case Euler::Explicit:
            // Explicit (Forward Euler)
            return this->EuroCallDiscreteQ_expl(M_1);
            break;
        case Euler::ImEx:
            // Implicit-explicit (Crank-Nicolson)
            return this->EuroCallDiscreteQ_imex(M_1);
            break;
        default:
            // Default to explicit method
            return this->EuroCallDiscreteQ_expl(M_1);
            break;
    }
}

std::vector<double> FiniteDifferencePricer::EuroCallDiscreteQ_expl(std::size_t M_1) {
    
    std::vector<double> res;
    
    // 1. Computational domain (part 2)
    // M: number of time intervals on the tau-axis
    // N: number of x intervals on the x-axis
    std::size_t N, M_2; double dtau_1, dtau_2, dx, alpha_2, x_compute, x_bar_compute;
    std::tie(N, M_2, dx, dtau_1, dtau_2, alpha_2, x_compute, x_bar_compute) = this->DomainParams_discrete_divs(M_1);
    
    // 3. Finite difference scheme
    // Fill x mesh
    std::vector<double> x_mesh;
    x_mesh = this->BuildMesh_discrete_divs(N, dx);
    
    std::vector<std::vector<double>> u_meshes = this->FiniteDifference_expl_discr(x_mesh, M_1, dtau_1, M_2, dtau_2, alpha_2);
    std::vector<double> u_mesh = u_meshes[0];
    std::vector<double> u_mesh_prev = u_meshes[1];
    std::transform(x_mesh.begin(), x_mesh.end(), x_mesh.begin(), std::bind(std::minus<double>(), std::placeholders::_1, std::log(1 - q_)));
    
    res.insert(res.end(), u_mesh[idx_compute]);
    
    double V_approx = std::exp(-a_ * x_mesh[idx_compute] - b_ * tau_final_) * u_mesh[idx_compute];
    
    res.insert(res.end(), V_approx);
    
    // 4. Pointwise convergence
    //std::vector<double> approximations = this->Approximate(x_mesh, u_mesh, dx);
    //res.insert(res.end(), approximations.cbegin(), approximations.cend());
    //double V_approx = approximations.front();
    
    
    // 5. RMS error
    //double error_RMS = this->RMS(x_mesh, u_mesh);
    //res.push_back(error_RMS);
    
    // 6. Greeks
    std::vector<double> Greeks = this->Greeks(idx_compute, x_mesh, u_mesh, u_mesh_prev, dtau_2, res[1]);
    res.insert(res.end(), Greeks.cbegin(), Greeks.cend());
    
    return res;
}

std::vector<double> FiniteDifferencePricer::EuroCallDiscreteQ_imex(std::size_t M_1) {
    
    std::vector<double> res;
    
    // 1. Computational domain (part 2)
    // M: number of time intervals on the tau-axis
    // N: number of x intervals on the x-axis
    std::size_t N, M_2; double dtau_1, dtau_2, dx, alpha_2, x_compute, x_bar_compute;
    std::tie(N, M_2, dx, dtau_1, dtau_2, alpha_2, x_compute, x_bar_compute) = this->DomainParams_discrete_divs(M_1);
    
    // 3. Finite difference scheme
    // Fill x mesh
    std::vector<double> x_mesh;
    x_mesh = this->BuildMesh_discrete_divs(N, dx);
    
    std::vector<std::vector<double>> u_meshes = this->FiniteDifference_discrete_imex(x_mesh, M_1, dtau_1, M_2, dtau_2, alpha_2);
    std::vector<double> u_mesh = u_meshes[0];
    std::vector<double> u_mesh_prev = u_meshes[1];
    
    std::transform(x_mesh.begin(), x_mesh.end(), x_mesh.begin(), std::bind(std::minus<double>(), std::placeholders::_1, std::log(1 - q_)));
    
    res.insert(res.end(), u_mesh[idx_compute]);
    
    
    
    double V_approx = std::exp(-a_ * x_compute - b_ * tau_final_) * u_mesh[idx_compute];
    
    res.insert(res.end(), V_approx);
    
    // 4. Pointwise convergence
    //std::vector<double> approximations = this->Approximate(x_mesh, u_mesh, dx);
    //res.insert(res.end(), approximations.cbegin(), approximations.cend());
    //double V_approx = approximations.front();
    
    
    // 5. RMS error
    //double error_RMS = this->RMS(x_mesh, u_mesh);
    //res.push_back(error_RMS);
    
    // 6. Greeks
    std::vector<double> Greeks = this->Greeks(idx_compute, x_mesh, u_mesh, u_mesh_prev, dtau_2, res[1]);
    res.insert(res.end(), Greeks.cbegin(), Greeks.cend());
    
    return res;
}

std::vector<double> FiniteDifferencePricer::DownAndOut(std::size_t M_1, const Euler& method){
    
    // 2. Boundary conditions
    boundary_tau_0_ = [=](double x)->double {
        // K * exp(ax) * max(1 - exp(x), 0)
        if (x > 0.) {
            return K_ * std::exp(a_ * x) * (std::exp(x) - 1.);
        } else {
            return 0.;
        }
    };
    
    boundary_x_l_ = [](double tau, double x)->double {
        return 0.;
    };
    
    boundary_x_r_ = [=](double tau, double x)->double {
    
        return K_ * std::exp(a_ * x + b_ * tau) * (std::exp(x - 2. * q_ * tau / (sigma_ * sigma_) ) - std::exp(-2. * r_ * tau / (sigma_ * sigma_)));
    };
    
    switch (method) {
        case Euler::Explicit:
            // Explicit (Forward Euler)
            return this->DownAndOut_expl(M_1);
            break;
        case Euler::Implicit:
            // Explicit (Forward Euler)
            return this->DownAndOut_impl(M_1);
            break;
        case Euler::ImEx:
            // Implicit-explicit (Crank-Nicolson)
            return this->DownAndOut_imex(M_1);
            break;
        default:
            // Default to explicit method
            return this->DownAndOut_expl(M_1);
            break;
    }
}


std::vector<double> FiniteDifferencePricer::DownAndOut_expl(std::size_t M) {
    
    std::vector<double> res;
    
    // 1. Computational domain (part 2)
    // M: number of time intervals on the tau-axis
    // N: number of x intervals on the x-axis
    std::size_t N, interval_i; double dtau, dx, alpha, x_compute;
    std::tie(N, dx, dtau, alpha, x_compute) = this->DomainParams_barrier(M);
    
    // 3. Finite difference scheme
    // Fill x mesh
    std::vector<double> x_mesh;
    std::tie(x_mesh, interval_i) = this->BuildMesh(N, dx);

    std::vector<std::vector<double>> u_meshes = this->FiniteDifference(alpha, x_mesh, M, dtau);
    std::vector<double> u_mesh = u_meshes[0];
    std::vector<double> u_mesh_prev = u_meshes[1];
    
    res.insert(res.end(), u_mesh[interval_i]);
    
    // 4. Pointwise convergence
    std::vector<double> approximations = this->Approximate(x_mesh, u_mesh, dx);
    //res.insert(res.end(), approximations.cbegin(), approximations.cend());
    double V_approx = approximations.front();
    res.insert(res.end(), V_approx);
    
    double error = std::abs(V_approx - 4.21963759) ;
    res.push_back(error);
    

    // 6. Greeks
    std::vector<double> Greeks = this->Greeks(interval_i, x_mesh, u_mesh, u_mesh_prev, dtau, V_approx);
    res.insert(res.end(), Greeks.cbegin(), Greeks.cend());
    
    return res;
}

std::vector<double> FiniteDifferencePricer::DownAndOut_impl(std::size_t M) {
    
    std::vector<double> res;
    
    // 1. Computational domain (part 2)
    // M: number of time intervals on the tau-axis
    // N: number of x intervals on the x-axis
    std::size_t N, interval_i; double dtau, dx, alpha, x_compute;
    std::tie(N, dx, dtau, alpha, x_compute) = this->DomainParams_barrier(M);
    
    // 3. Finite difference scheme
    // Fill x mesh
    std::vector<double> x_mesh;
    std::tie(x_mesh, interval_i) = this->BuildMesh(N, dx);

    std::vector<std::vector<double>> u_meshes = this->FiniteDifference_impl(alpha, x_mesh, M, dtau);
    std::vector<double> u_mesh = u_meshes[0];
    std::vector<double> u_mesh_prev = u_meshes[1];
    
    res.insert(res.end(), u_mesh[interval_i]);
    
    // 4. Pointwise convergence
    std::vector<double> approximations = this->Approximate(x_mesh, u_mesh, dx);
    //res.insert(res.end(), approximations.cbegin(), approximations.cend());
    double V_approx = approximations.front();
    res.insert(res.end(), V_approx);
    
    double error = std::abs(V_approx - 4.21963759) ;
    res.push_back(error);
    
    // 6. Greeks
    std::vector<double> Greeks = this->Greeks(interval_i, x_mesh, u_mesh, u_mesh_prev, dtau, V_approx);
    res.insert(res.end(), Greeks.cbegin(), Greeks.cend());
    
    return res;
}

std::vector<double> FiniteDifferencePricer::DownAndOut_imex(std::size_t M) {
    
    std::vector<double> res;
    
    // 1. Computational domain (part 2)
    // M: number of time intervals on the tau-axis
    // N: number of x intervals on the x-axis
    std::size_t N, interval_i; double dtau, dx, alpha, x_compute;
    std::tie(N, dx, dtau, alpha, x_compute) = this->DomainParams_barrier(M);
    
    // 3. Finite difference scheme
    // Fill x mesh
    std::vector<double> x_mesh;
    std::tie(x_mesh, interval_i) = this->BuildMesh(N, dx);

    std::vector<std::vector<double>> u_meshes = this->FiniteDifference_imex(alpha, x_mesh, M, dtau);
    std::vector<double> u_mesh = u_meshes[0];
    std::vector<double> u_mesh_prev = u_meshes[1];
    
    res.insert(res.end(), u_mesh[interval_i]);
    
    // 4. Pointwise convergence
    std::vector<double> approximations = this->Approximate(x_mesh, u_mesh, dx);
    //res.insert(res.end(), approximations.cbegin(), approximations.cend());
    double V_approx = approximations.front();
    res.insert(res.end(), V_approx);
    
    double error = std::abs(V_approx - 4.21963759) ;
    res.push_back(error);
    
    // 6. Greeks
    std::vector<double> Greeks = this->Greeks(interval_i, x_mesh, u_mesh, u_mesh_prev, dtau, V_approx);
    res.insert(res.end(), Greeks.cbegin(), Greeks.cend());
    
    return res;
}

