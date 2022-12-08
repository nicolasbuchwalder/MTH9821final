//
//  LinearSolver.hpp
//  MTH9815_HW12_FD
//
//  Created by Nicolas Buchwalder on 28.11.22.
//

#ifndef LinearSolver_hpp
#define LinearSolver_hpp

#include "EigenCommonHeader.h"

class LinearSolver {
public:
    vec ForwardSub(const mat& L, const vec& b) const;
    vec BackwardSub(const mat& U, const vec& b) const;
    vec LUSolve(const mat& A, const vec& b) const;
    vec CholeskySolve(const mat& A, const vec& b) const;
    
    vec DiagonalSolve_vec(const mat& D, const vec& b) const;
    mat DiagonalSolve(const mat& D, const mat& B) const;
    
    mat ForwardSub_mat(const mat& L, const mat& B) const;
};

#endif /* LinearSolver_hpp */
