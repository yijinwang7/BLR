// Copyright (c) 2025 Yijin Wang
// This software is released under the MIT License.
// https://opensource.org/licenses/MIT

#ifndef BLOCKLLLHYBKZ_H
#define BLOCKLLLHYBKZ_H

#include <Eigen/Dense>
#include <vector>
#include "randomPermutation.h"
#include "lllmatrix.h"
#include <fplll.h>
#include "utils.h"
#include "BLR.h"  // Provides definition of Config and LLLMatrix

namespace BLR {
    LLLMatrix blockLLLHy_BKZ(LLLMatrix B, int p, int fre, int b_size, double stopC, const Config &config);
}

#endif // BLOCKLLLHYBKZ_H
