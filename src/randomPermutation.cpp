// Copyright (c) 2025 Yijin Wang
// This software is released under the MIT License.
// https://opensource.org/licenses/MIT

#include "randomPermutation.h"

//LLLMatrix randomPermutation(const LLLMatrix &B)
//{
//    int m = B.cols(); // Assuming B is a matrix, get the number of columns
//
//    // Check if we need to apply a random permutation after a full sweep
//
//    // Create a vector for indices
//    std::vector<int> perm(m);
//    for (int i = 0; i < m; ++i)
//    {
//        perm[i] = i;
//    }
//
//    // Shuffle the indices
//    std::srand(static_cast<unsigned int>(std::time(0))); // Seed with current time
//    std::random_shuffle(perm.begin(), perm.end());
//
//    // Apply the permutation to the matrix
//    LLLMatrix B_permuted = B;
//    for (int i = 0; i < m; ++i)
//    {
//        B_permuted.col(i) = B.col(perm[i]);
//    }
//
//    // Copy the permuted matrix back to B
//    // B = B_permuted;
//
//    return B_permuted;
//}

void randomPermutation2(LLLMatrix &B,std::mt19937 &g)
{
    int m = B.cols(); // Assuming B is a matrix, get the number of columns
    //int rows = B.rows();

    // Check if we need to apply a random permutation after a full sweep

    // Create a vector for indices
    std::vector<int> perm(m);
    for (int i = 0; i < m; ++i)
    {
        perm[i] = i;
    }

    std::shuffle(perm.begin(), perm.end(), g);

    LLLMatrix B_temp = B;

    // Apply the permutation to the matrix
    for (int i = 0; i < m; i++)
    {
        B.col(i) = B_temp.col(perm[i]);
    }
}

//void randomPermutation3(LLLMatrix &B, fplll::ZZ_mat<mpz_t> &U) {
//    int m = B.cols(); // Assuming B is a matrix, get the number of columns
//    //int rows = B.rows(); // Get the number of rows
//
//    // Initialize U as an m x m identity matrix
//    U.resize(m, m); // Ensure U is the correct size
//    for (int i = 0; i < m; i++) {
//        for (int j = 0; j < m; j++) {
//            U[i][j] = (i == j) ? 1 : 0;
//        }
//    }
//
//    // Create a vector for indices
//    std::vector<int> perm(m);
//    for (int i = 0; i < m; ++i) {
//        perm[i] = i;
//    }
//
//    // Shuffle the indices
//    std::srand(static_cast<unsigned int>(std::time(0))); // Seed with current time
//    std::random_shuffle(perm.begin(), perm.end());
//
//    // Temporary copy to hold original columns of B
//    LLLMatrix B_temp = B;
//    LLLMatrix U_temp = zzmat2lllmatrix(U);
//    LLLMatrix U_final = U_temp;
//
//    // Apply the permutation to the matrix B and U
//    for (int i = 0; i < m; i++) {
//        B.col(i) = B_temp.col(perm[i]); // Assuming col() gets and setCol() sets the column in B
//        U_final.col(i) = U_temp.col(perm[i]);
////        for (int j = 0; j < m; j++) {
////            U[j][i] = (perm[i] == j) ? 1 : 0; // Set column i of U according to perm
////        }
//
//    }
//    U = lllmatrix2zzmat(U_final);
//}
