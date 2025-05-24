// Copyright (c) 2025 Yijin Wang
// This software is released under the MIT License.
// https://opensource.org/licenses/MIT

#ifndef LINEARSYSTEMSOLVER_RANDOMPERMUTATION_H
#define LINEARSYSTEMSOLVER_RANDOMPERMUTATION_H
#include <Eigen/Dense>
#include <random>
#include <algorithm> // for std::shuffle
#include <vector>
#include "lllmatrix.h"
#include <fplll.h>
#include "utils.h"

//LLLMatrix randomPermutation(const LLLMatrix &B);
void randomPermutation2(LLLMatrix &B, std::mt19937 &g);
//void randomPermutation3(LLLMatrix &B, fplll::ZZ_mat<mpz_t> &U);
#endif // LINEARSYSTEMSOLVER_RANDOMPERMUTATION_H

