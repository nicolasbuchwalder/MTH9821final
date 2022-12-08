//
//  FiniteDifferencePricer.hpp
//  MTH9815_HW12_FD
//
//  Created by Nicolas Buchwalder on 28.11.22.
//

#ifndef FiniteDifferencePricer_hpp
#define FiniteDifferencePricer_hpp

#include "Option.hpp"
#include <functional>
#include <vector>
#include <tuple>
#include <array>
#include "EigenCommonHeader.h"

enum Euler {
    Explicit,
    Implicit,
    ImEx,
    CN_LU,
    CN_SOR,
};

class FiniteDifferencePricer {
private:
    // Option data
    double S0_;
    double K_;
    double T_;
    double sigma_;
    double r_;
    double q_;
    double t_div_;
    double _B;
    
    // Heat equation transformation coefficients
    double a_;
    double b_;
    
    // Finite difference parameters
    double tau_div_;
    double tau_final_;
    double x_l_;
    double x_r_;
    double x_l_new_;
    double x_r_new_;
    
    std::size_t idx_compute;
    
//    double dt_;
//    double dx_;
//    std::size_t N_;
//    double alpha_;
    
    // Boundary conditions
    std::function<double (double)> boundary_tau_0_;
    std::function<double (double, double)> boundary_x_l_;
    std::function<double (double, double)> boundary_x_r_;
    
    
    // 1. Computational domain
    std::tuple<std::size_t, double, double, double> DomainParams(std::size_t M) const;
    std::tuple<std::size_t, double, double, double, double, double, double, double> DomainParams_discrete_divs(std::size_t M_1);
    std::tuple<size_t, double, double, double, double> DomainParams_barrier(std::size_t M);
    
    std::pair<std::vector<double>, std::size_t> BuildMesh(std::size_t N, double dx) const;
    
    std::vector<double> BuildMesh_discrete_divs(std::size_t N, double dx) const;
    
    // 3. Finite difference scheme
    std::vector<std::vector<double>> FiniteDifference(double alpha, const std::vector<double>& x_mesh, std::size_t M, double dtau) const;
    
    std::vector<std::vector<double>> FiniteDifference_AmeriPut(double alpha, const std::vector<double>& x_mesh, std::size_t M, double dtau) const;
    
    std::vector<std::vector<double>> FiniteDifference_impl(double alpha, const std::vector<double>& x_mesh, std::size_t M, double dtau) const;
    
    std::vector<std::vector<double>> FiniteDifference_imex(double alpha, const std::vector<double>& x_mesh, std::size_t M, double dtau) const;
    
    std::vector<std::vector<double>> FiniteDifference_imex_AmeriPut(double alpha, const std::vector<double>& x_mesh, std::size_t M, double dtau) const;
    
    std::vector<std::vector<double>> FiniteDifference_expl_discr(const std::vector<double>& x_mesh, std::size_t M_1, double dtau_1, std::size_t M_2, double dtau_2, double alpha_2) const;
    
    std::vector<std::vector<double>> FiniteDifference_discrete_imex(const std::vector<double>& x_mesh, std::size_t M_1, double dtau_1, std::size_t M_2, double dtau_2, double alpha_2) const;
    
    std::vector<double> EuroCallDiscreteQ_expl(std::size_t M_1);
    std::vector<double> EuroCallDiscreteQ_imex(std::size_t M_1);
    
    std::vector<double> DownAndOut_expl(std::size_t M);
    std::vector<double> DownAndOut_impl(std::size_t M);
    std::vector<double> DownAndOut_imex(std::size_t M);
    
    // Advance time (modify u-mesh in place)
    // Forward Euler
    void EuroPut_advance(double tau, double xl, double xr, double alpha, const std::vector<double>& x_mesh, std::vector<double>& u_mesh) const;
    void AmeriPut_advance(double tau, double alpha, const std::vector<double>& x_mesh, std::vector<double>& u_mesh) const;
    
    void EuroCallDiscrete_advance(double tau1, double tau2);
    
    // Backward Euler
    void EuroPut_advance_implicit(double tau, double alpha, const std::vector<double>& x_mesh, std::vector<double>& u_mesh, const mat& L, const mat& U) const;
    
    // Crank-Nicolson
    void EuroPut_advance_imex(double tau, double xl, double xr, double alpha, const std::vector<double>& x_mesh, std::vector<double>& u_mesh, const mat& L, const mat& U, const mat& b_multiplier) const;
    void AmeriPut_advance_imex(double tau, double alpha, const std::vector<double>& x_mesh, std::vector<double>& u_mesh, const mat& A, const mat& b_multiplier) const;
    
    // 4. Pointwise convergence
    std::vector<double> Approximate(const std::vector<double>& x_mesh, const std::vector<double>& u_mesh, double dx) const;
    
    // 5. RMS error
    double RMS(const std::vector<double>& x_mesh, const std::vector<double>& u_mesh) const;
    
    // 6. Greeks
    std::vector<double> Greeks(std::size_t interval_i, const std::vector<double>& x_mesh, const std::vector<double>& u_mesh, std::vector<double>& u_mesh_prev, double dtau, double V_approx) const;
    
    // 7. Variance reduction
    double VarianceReduction_AmeriPut(double V_approx, std::size_t M) const;
    
    double VarianceReduction_imex_AmeriPut(double V_approx, std::size_t M) const;
    
public:
    FiniteDifferencePricer(double S0, double K, double T, double sigma, double r, double q, double t_div, double B);
    ~FiniteDifferencePricer() = default;
    void PrintVector(const std::vector<double>& vec) const;
    
    // Finite difference hyperparameters
    static double alpha_temp_;
    static std::size_t M_init_;
    
    std::vector<double> EuroPut(std::size_t M, const Euler& method = Euler::Explicit);
    // Explicit method (Forward Euler)
    std::vector<double> EuroPut_expl(std::size_t M) const;
    // Implicit method (Backward Euler)
    std::vector<double> EuroPut_impl(std::size_t M) const;
    // ImEx method (Crank-Nicolson)
    std::vector<double> EuroPut_imex(std::size_t M) const;
    
    std::vector<double> AmericanPut(std::size_t M, const Euler& method = Euler::Explicit);
    std::vector<double> AmericanPut_expl(std::size_t M) const;
    std::vector<double> AmericanPut_imex(std::size_t M) const;
    
    std::vector<double> EuroCallDiscreteQ(std::size_t M_1, const Euler& method);
    
    std::vector<double> DownAndOut(std::size_t M_1, const Euler& method);
    
    std::array<std::vector<double>, 2> AmericanPut_EarlyExDomain(std::size_t M);
};

#endif /* FiniteDifferencePricer_hpp */
