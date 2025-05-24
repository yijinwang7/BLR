// Copyright (c) 2025 Yijin Wang
// This software is released under the MIT License.
// https://opensource.org/licenses/MIT

#ifndef BLOCKLLLHYSTANDARD_H
#define BLOCKLLLHYSTANDARD_H

#include <Eigen/Dense>
#include <vector>
#include "randomPermutation.h"
#include "lllmatrix.h"
#include <fplll.h>
#include "utils.h"
#include "BLR.h"  // Provides definition of Config and LLLMatrix

namespace BLR {

    // Declaration for the standard version of blockLLLHy.
    // Accepts a lattice matrix, number of blocks (p), frequency (fre),
    // stopping criteria (stopC), and the configuration.
    LLLMatrix blockLLLHy_standard(LLLMatrix B, int p, int fre, double stopC, const Config &config);

} // namespace BLR

#endif // BLOCKLLLHYSTANDARD_H
