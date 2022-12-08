//
//  FiniteDifference.hpp
//  FD
//
//  Created by Nicolas Buchwalder on 07.12.22.
//

#ifndef FiniteDifference_hpp
#define FiniteDifference_hpp

#include "EigenCommonHeader.h"
#include "Option.hpp"

#include <functional>
#include <vector>
#include <tuple>
#include <array>

enum Euler {
    Expl,
    Impl,
    CN_LU,
    CN_SOR,
};

class FiniteDifference {

private:
    
    // flags
    bool _params_called;    // flag for update_params being called
    bool _boundaries_called;// flag for set_boundaries being called
    
    OptionPayoff _payoff;   // payoff type
    OptionExercise _ex;     // exercise type
    OptionType _type;       // type of option
    double _S;              // Spot price
    double _K;              // Strike price
    double _T;              // Maturity
    double _sigma;          // Volatility
    double _r;              // Const interest rate
    double _q;              // Compound dividend rate
    
    DivsTuple _divs;        // dividend
    std::vector<double> _add_params;    // additional parameters of option
    
    
    // domain coefficients
    std::vector<double> _Ms;
    std::vector<double> _Ns;
    
    std::vector<double> _x_ls;
    std::vector<double> _x_rs;
    std::vector<double> _taus;
    
    std::vector<double> _alphas;
    std::vector<double> _dxs;
    std::vector<double> _dtaus;
    
    // heat coefficients
    double _a;
    double _b;
    double _tau_final_;
    
    
    // meshes
    std::vector<double> _x_mesh;
    std::vector<std::vector<double>> _u_mesh;
    
    
    // boundary conditions
    std::function<double (double)> _boundary_tau_0_;
    std::function<double (double, double)> _boundary_x_l_;
    std::function<double (double, double)> _boundary_x_r_;
    
    
    void build_mesh();
    
    // functions to compute advance one step
    void advance_expl(double alpha, double tau, double xl, double xr);
    void advance_impl(double alpha, double tau, double xl, double xr, const mat& L, const mat& U);
    void advance_cn_lu();
    void advance_cn_sor();
    
    // functions to compute the finite difference over subdomain
    void compute_sub_domain_expl();
    void compute_sub_domain_impl();
    void compute_sub_domain_cn_lu();
    void compute_sub_domain_cn_sor();
    
    // functions to compute all the
    void price_expl();
    void price_impl();
    void price_cn_lu();
    void price_cn_sor();
    
    
public:
    
    FiniteDifference() = default;
    
    // set all parameters
    void update_params(Option opt, std::size_t M_1, double alpha_1_temp);
    
    // price option
    std::vector<double> price_option(const Euler& scheme, bool include_greeks);
    
    // print the approximations
    void show_grid();
    

};

#endif /* FiniteDifference_hpp */
