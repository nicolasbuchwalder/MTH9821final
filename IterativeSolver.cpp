//
//  IterativeSolver.cpp
//  MTH9815_HW12_FD
//
//  Created by Nicolas Buchwalder on 28.11.22.
//

#include "IterativeSolver.hpp"
#include "LinearSolver.hpp"

IterativeSolver::IterativeSolver(const mat& A, const vec& b) : A_(A), b_(b) {
    x0_ = vec::Zero(b.size());
}

IterativeSolver::IterativeSolver(const mat& A, const vec& b, const vec& x0) : A_(A), b_(b), x0_(x0) {}

void IterativeSolver::x0(const vec& init) {
    x0_ = init;
}

std::tuple<mat, mat, mat> IterativeSolver::LDU(const mat& A) const {
    mat D = A.diagonal().asDiagonal();
    mat L = A.triangularView<Eigen::StrictlyLower>();
    mat U = A.triangularView<Eigen::StrictlyUpper>();
    
    return std::make_tuple(L, D, U);
}


std::tuple<mat, vec> IterativeSolver::Jacobi_RC(const mat& A, const vec& b) const {
    
    mat L, D, U;
    std::tie(L, D, U) = this->LDU(A);
    
    LinearSolver linear;
    vec c(linear.DiagonalSolve_vec(D, b));
    mat R(-linear.DiagonalSolve(D, L + U));
    
    return std::make_tuple(R, c);
}

std::tuple<mat, vec> IterativeSolver::GaussSeidel_RC(const mat& A, const vec& b) const {
    
    mat L, D, U;
    std::tie(L, D, U) = this->LDU(A);
    
    LinearSolver linear;
    vec c(linear.ForwardSub(D + L, b));
    mat R(-linear.ForwardSub_mat(D + L, U));
    
    return std::make_tuple(R, c);
}

std::tuple<mat, vec> IterativeSolver::SOR_RC(const mat& A, const vec& b, double omega) const {
    
    mat L, D, U;
    std::tie(L, D, U) = this->LDU(A);
    
    LinearSolver linear;
    vec c(linear.ForwardSub(D + omega * L, omega * b));
    mat R(linear.ForwardSub_mat(D + omega * L, (1. - omega) * D - omega * U));
    
    return std::make_tuple(R, c);
}



std::tuple<vec, unsigned> IterativeSolver::Recursion(const vec& x0, const mat& R, const vec& c, const mat& A, const vec& b, const StoppingCriterion& criterion, const double tolerance) const {
    
    unsigned iter_N = 0;
    vec xold(x0);
    
    switch (criterion) {
        case StoppingCriterion::consecutive:
        {
            // Consecutive approximation criterion
            vec diff(vec::Ones(xold.size()));
            vec xnew(xold.size());
            while (diff.norm() > tolerance) {
                xnew = R * xold + c;
                diff = xnew - xold;
                xold = xnew;
                iter_N++;
            }
            break;
        }
        case StoppingCriterion::residual:
        {
            // Residual-based criterion
            vec r(b - A * xold);
            vec xnew(xold.size());
            double residual_criterion = tolerance * r.norm();
            while (r.norm() > residual_criterion) {
                xnew = R * xold + c;
                r = b - A * xnew;
                xold = xnew;
                iter_N++;
            }
            break;
        }
            //        case StoppingCriterion::combined:
            //        {
            //            // Combined criteria
            //            vec r(b - A * xold);
            //            vec diff(xold.Identity());
            //            vec xnew(xold.size());
            //            double residual_criterion = tolerance * r.norm();
            //            while (r.norm() > residual_criterion) {
            //                xnew = R * xold + c;
            //                r = b - A * xold;
            //                xold = xnew;
            //            }
            //            break;
            //        }
        default:
        {
            // Do nothing
        }
            break;
    }
    
    return std::make_tuple(xold, iter_N);
}

std::tuple<vec, unsigned> IterativeSolver::RecursionProjected_lower(const vec& x0, const mat& R, const vec& c, const mat& A, const vec& b, const StoppingCriterion& criterion, const double tolerance, const vec& bound) const {
    
    unsigned iter_N = 0;
    vec xold(x0);
    
    switch (criterion) {
        case StoppingCriterion::consecutive:
        {
            // Consecutive approximation criterion
            vec diff(vec::Ones(xold.size()));
            vec xnew(xold.size());
            while (diff.norm() > tolerance) {
                xnew = R * xold + c;
                for (std::size_t i = 0; i < xnew.size(); i++) {
                    xnew(i) = std::max(xnew(i), bound(i));
                }
                diff = xnew - xold;
                xold = xnew;
                iter_N++;
            }
            break;
        }
        case StoppingCriterion::residual:
        {
            // Residual-based criterion
            vec r(b - A * xold);
            vec xnew(xold.size());
            double residual_criterion = tolerance * r.norm();
            while (r.norm() > residual_criterion) {
                xnew = R * xold + c;
                for (std::size_t i = 0; i < xnew.size(); i++) {
                    xnew(i) = std::max(xnew(i), bound(i));
                }
                r = b - A * xnew;
                xold = xnew;
                iter_N++;
            }
            break;
        }
            //        case StoppingCriterion::combined:
            //        {
            //            // Combined criteria
            //            vec r(b - A * xold);
            //            vec diff(xold.Identity());
            //            vec xnew(xold.size());
            //            double residual_criterion = tolerance * r.norm();
            //            while (r.norm() > residual_criterion) {
            //                xnew = R * xold + c;
            //                r = b - A * xold;
            //                xold = xnew;
            //            }
            //            break;
            //        }
        default:
        {
            // Do nothing
        }
            break;
    }
    
    return std::make_tuple(xold, iter_N);
}



std::tuple<vec, unsigned> IterativeSolver::Jacobi(const StoppingCriterion& criterion, double tolerance) const {
    
    mat R;
    vec c;
    std::tie(R, c) = Jacobi_RC(A_, b_);
    
    return this->Recursion(x0_, R, c, A_, b_, criterion, tolerance);
}

std::tuple<vec, unsigned> IterativeSolver::GaussSeidel(const StoppingCriterion& criterion, double tolerance) const {
    
    mat R;
    vec c;
    std::tie(R, c) = GaussSeidel_RC(A_, b_);
    
    return this->Recursion(x0_, R, c, A_, b_, criterion, tolerance);
}

std::tuple<vec, unsigned> IterativeSolver::SOR(double omega, const StoppingCriterion& criterion, double tolerance) const {
    
    mat R;
    vec c;
    std::tie(R, c) = SOR_RC(A_, b_, omega);
    
//    std::cout << "RC" << std::endl;
//    std::cout << R << std::endl;
//    std::cout << c << std::endl;
    
    return this->Recursion(x0_, R, c, A_, b_, criterion, tolerance);
}

std::tuple<vec, unsigned> IterativeSolver::SORProjected_lower(double omega, const StoppingCriterion& criterion, double tolerance, const vec& bound) const {
    
    mat R;
    vec c;
    std::tie(R, c) = SOR_RC(A_, b_, omega);
    
//    std::cout << "RC" << std::endl;
//    std::cout << R << std::endl;
//    std::cout << c << std::endl;
    
    return this->RecursionProjected_lower(x0_, R, c, A_, b_, criterion, tolerance, bound);
}

