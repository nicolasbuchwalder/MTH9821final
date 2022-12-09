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
#include <iostream>

enum Scheme {
    eul_expl,
    eul_impl,
    cn_lu,
    cn_sor,
};

class FiniteDifference {

private:
    
    Option _option;         // option
    OptionPayoff _payoff;   // payoff type
    OptionExercise _ex;     // exercise type
    OptionType _type;       // type of option
    double _S;              // Spot price
    double _K;              // Strike price
    double _T;              // Maturity
    double _sigma;          // Volatility
    double _r;              // Const interest rate
    double _q;              // Compound dividend rate
    
    std::size_t _num_divs;  // number of discrete dividends
    DivsTuple _divs;        // discrete dividends
    std::vector<double> _add_params;    // additional parameters of option
    
    
    // domain coefficients
    std::vector<std::size_t> _Ms;
    std::size_t _N;
    double _alpha_temp;
    
    double _x_l;
    double _x_r;
    std::vector<double> _taus;
    
    std::vector<double> _alphas;
    std::vector<double> _dxs;
    std::vector<double> _dtaus;
    
    // heat coefficients
    double _a;
    double _b;
    double _tau_final;
    
    // list of dividends
    std::vector<double> _tau_divs;
    std::vector<double> _q_divs;
    
    // meshes
    std::vector<double> _x_mesh;
    std::vector<std::vector<double>> _u_mesh;
    
    // terminal values
    double _x_compute;
    std::pair<std::size_t, std::size_t> _x_compute_idx;
    
    
    // boundary conditions
    std::function<double (double)> _boundary_tau_0;
    std::function<double (double, double)> _boundary_x_l;
    std::function<double (double, double)> _boundary_x_r;
    
    // domain builders
    void build_domain();
    // sub functions
    void set_domain();  // setting the
    void set_discretisation();  // setting the N, alpha, dx, dtaus etc
    void set_boundaries(); // setting the boundary functions in function of type of option
    void build_mesh();  // constructing the mesh
    
    
    // functions to compute advance one step
    void advance_expl(double alpha, double tau, double xl, double xr);
    void advance_impl(double alpha, double tau, double xl, double xr, const mat& L, const mat& U);
    void advance_cn_lu(double alpha, double tau, double xl, double xr, const mat& L, const mat& U, const mat& b_multiplier);
    void advance_cn_sor();
    
    // functions to compute the finite difference over subdomain (after each discrete divs)
    void compute_sub_domain_expl(std::size_t sub, double start_tau);
    void compute_sub_domain_impl(std::size_t sub, double start_tau);
    void compute_sub_domain_cn_lu(std::size_t sub, double start_tau);
    void compute_sub_domain_cn_sor(std::size_t sub, double start_tau);
    
    // function to convert from heat to real world
    double convert_to_v(double x, double tau, double u) const;
    
    double approximate();
    
    // computing the greeks
    std::vector<double> greeks();
    
    
public:
    
    FiniteDifference() = default;
    
    // set all parameters
    void set_params(Option opt, std::size_t M_1, double alpha_temp);
    
    // price option
    std::vector<double> price_option(const Scheme& scheme, bool include_greeks);
    
    // show the domain parameters
    void show_domain_params();
    
    // print the approximations
    void show_grid(bool convert);
    
    
    // print functions
    template < typename T >
    void print(const std::vector<T>& vec) const;
    template < typename T >
    void print(const std::vector<std::vector<T>>& mat) const;
    

};

#endif /* FiniteDifference_hpp */
