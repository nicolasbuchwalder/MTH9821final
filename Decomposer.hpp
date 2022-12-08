//
//  Decomposer.hpp
//  MTH9815_HW12_FD
//
//  Created by Nicolas Buchwalder on 28.11.22.
//

#ifndef Decomposer_hpp
#define Decomposer_hpp

#include <tuple>

#include "EigenCommonHeader.h"

class Decomposer {
   public:
    std::tuple<mat, mat> lu_no_pivoting(const mat& A) const;
    std::tuple<mat, mat, mat> lu_row_pivoting(const mat& A) const;
    mat cholesky(const mat& A) const;
};

#endif /* Decomposer_hpp */
