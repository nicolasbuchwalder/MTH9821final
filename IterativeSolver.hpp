//
//  IterativeSolver.hpp
//  MTH9815_HW12_FD
//
//  Created by Nicolas Buchwalder on 28.11.22.
//

#ifndef IterativeSolver_hpp
#define IterativeSolver_hpp

#include "EigenCommonHeader.h"
#include <tuple>
#include <functional>

enum StoppingCriterion {consecutive, residual, combined};

class IterativeSolver {
private:
    mat A_;
    vec b_;
    vec x0_;
    
    std::tuple<mat, mat, mat> LDU(const mat& A) const;
    
    std::tuple<mat, vec> Jacobi_RC(const mat& A, const vec& b) const;
    std::tuple<mat, vec> GaussSeidel_RC(const mat& A, const vec& b) const;
    std::tuple<mat, vec> SOR_RC(const mat& A, const vec& b, double omega) const;
    
    std::tuple<vec, unsigned> Recursion(const vec& x0, const mat& R, const vec& c, const mat& A, const vec& b, const StoppingCriterion& criterion, const double tolerance) const;
    
    std::tuple<vec, unsigned> RecursionProjected_lower(const vec& x0, const mat& R, const vec& c, const mat& A, const vec& b, const StoppingCriterion& criterion, const double tolerance, const vec& bound) const;
    std::tuple<vec, unsigned> RecursionProjected_lowerelementwise(const vec &x0, const mat &R, const vec &c, const mat &A, const vec &b, const StoppingCriterion &criterion, const double tolerance, const vec &bound, double alpha) const;

public:
    IterativeSolver(const mat& A, const vec& b);
    IterativeSolver(const mat& A, const vec& b, const vec& x0);
    
    void x0(const vec& init);
    
    std::tuple<vec, unsigned> Jacobi(const StoppingCriterion& criterion, double tolerance) const;
    std::tuple<vec, unsigned> GaussSeidel(const StoppingCriterion& criterion, double tolerance) const;
    std::tuple<vec, unsigned> SOR(double omega, const StoppingCriterion& criterion, double tolerance) const;
    
    std::tuple<vec, unsigned> SORProjected_lower(double omega, const StoppingCriterion& criterion, double tolerance, const vec& bound) const;
    std::tuple<vec, unsigned> SORProjected_lowerelementwise(double omega, const StoppingCriterion& criterion, double tolerance, const vec& bound, double alpha) const;
};

#endif /* IterativeSolver_hpp */
