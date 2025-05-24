// Copyright (c) 2025 Yijin Wang
// This software is released under the MIT License.
// https://opensource.org/licenses/MIT

#include "BLR.h"
#include <fplll.h>
#include <Eigen/Dense>
#include <vector>
#include <omp.h>
#include <fstream>
#include <iostream>
#include <functional>
#include <chrono>
#include <random>
#include <algorithm>
#include <gmp.h>
#include "utils.h"


int main() {
    // Create matrix B.
    int d = 20;
    fplll::ZZ_mat<mpz_t> matrix(d, d);
    matrix.gen_trg(0.9);
    auto B22 = zzmat2lllmatrix(matrix);
    auto fplllA = lllmatrix2zzmat(B22);
    LLLMatrix B2_1 = B22.transpose();
    LLLMatrix B = B2_1;

    // Create a configuration object with default values.
    BLR::Config config;

    // Customize the configuration as needed.
    config.reductionMethod = BLR::ReductionMethod::LLL;  // or HLLL, BKZ
    config.lllMethod = BLR::LLLMethod::Fast;              // choose Fast, Wrapper, etc.
    config.numBlocks = 5;                                 // e.g., split into 5 blocks
    config.customFre = 1;                                 // leave at 0 to use default behavior
    config.stopCriteria = 0.999;                          // adjust if needed
    config.ftType = BLR::FTType::Double;
    config.precision = 0;                               // for MPFR precision if used
    //config.b_size = 10;

    // Now call the unified reduction function.
    LLLMatrix reducedB = BLR::blockBLRReduction(B, config);
    // Test if the matrix reduced by blockLLLHy is LLL reduced
    fplll::ZZ_mat<mpz_t> fplllB1 = lllmatrix2zzmat(reducedB.transpose());
    fplll::ZZ_mat<mpz_t> arg_u(fplllB1.get_rows(), fplllB1.get_rows()); // Transformation matrix
    fplll::ZZ_mat<mpz_t> arg_uinv_t(fplllB1.get_rows(),
                                    fplllB1.get_rows()); // Inverse transpose of the transformation matrix


    fplll::MatGSO<fplll::Z_NR<mpz_t>, fplll::FP_NR<mpfr_t>> M(fplllB1, arg_u, arg_uinv_t, 0);
    // one on success
    int status = fplll::is_lll_reduced(M, fplll::LLL_DEF_DELTA, fplll::LLL_DEF_ETA);

    std::cout<< "one on success : " << status << std::endl;

    return 0;
}
