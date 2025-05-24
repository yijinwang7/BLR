// Copyright (c) 2025 Yijin Wang
// This software is released under the MIT License.
// https://opensource.org/licenses/MIT

#include "utils.h"
#include <iostream>

LLLMatrix zzmat2lllmatrix(const fplll::ZZ_mat<mpz_t> &B)
{

    const auto cols = B.get_cols();
    const auto rows = B.get_rows();

    LLLMatrix R(rows, cols);
    for (auto i = 0; i < rows; i++)
    {
        for (auto j = 0; j < cols; j++)
        {
            // R(i, j) = mpz_get_si((B[i][j]).get_data());
            R(i, j) = (B[i][j]).get_data();
        }
    }

    return R;
}

fplll::ZZ_mat<mpz_t> lllmatrix2zzmat(const LLLMatrix &B)
{
    const auto cols = B.cols();
    const auto rows = B.rows();

    fplll::ZZ_mat<mpz_t> R(rows, cols);
    for (auto i = 0; i < rows; i++)
    {
        for (auto j = 0; j < cols; j++)
        {
            mpz_t v;
            B(i, j).get_mpz(v);
            // mpz_init(v);
            // mpz_set_si(v, B(i, j));
            // R[i][j] = fplll::Z_NR<mpz_t>(v);
            R[i][j] = fplll::Z_NR<mpz_t>(v);
            mpz_clear(v);
        }
    }

    return R;
}

///**
// * Multiplies two matrices of type fplll::ZZ_mat<mpz_t>.
// *
// * @param u1 The first matrix.
// * @param u2 The second matrix.
// * @return The result of multiplying u1 by u2.
// */
//fplll::ZZ_mat<mpz_t> multiply_matrices(const fplll::ZZ_mat<mpz_t>& u1, const fplll::ZZ_mat<mpz_t>& u2)
//{
//    int m = u1.get_rows(); // Assuming square matrices for simplicity
//    fplll::ZZ_mat<mpz_t> result(m, m);
//
//    for (int i = 0; i < m; i++) {
//        for (int j = 0; j < m; j++) {
//            mpz_t sum;
//            mpz_init(sum); // Initialize sum to zero
//            for (int k = 0; k < m; k++) {
//                mpz_t product;
//                mpz_init(product);
//                mpz_mul(product, u1[i][k].get_data(), u2[k][j].get_data()); // Exact multiplication of mpz_t types
//                mpz_add(sum, sum, product); // Add the product to the sum for this element
//                mpz_clear(product);
//            }
//            result[i][j] = fplll::Z_NR<mpz_t>(sum); // Set the computed sum to the result matrix
//            mpz_clear(sum);
//        }
//    }
//
//    return result;
//}

/**
 * Multiplies two matrices of type fplll::ZZ_mat<mpz_t>, where the first is m x n and the second is n x n.
 *
 * @param u1 The first matrix, size m x n.
 * @param u2 The second matrix, size n x n.
 * @return The result of multiplying u1 by u2, size m x n.
 */
fplll::ZZ_mat<mpz_t> multiply_matrices(const fplll::ZZ_mat<mpz_t>& u1, const fplll::ZZ_mat<mpz_t>& u2)
{
    int m = u1.get_rows(); // Rows of the first matrix.
    int n = u2.get_cols(); // Columns of the second matrix, also the shared dimension.
    fplll::ZZ_mat<mpz_t> result(m, n); // Resultant matrix should be m x n.

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) { // Iterate up to n for the second matrix's columns.
            mpz_t sum;
            mpz_init(sum); // Initialize sum to zero.
            for (int k = 0; k < n; k++) { // Iterate over the shared dimension, which is now n.
                mpz_t product;
                mpz_init(product);
                mpz_mul(product, u1[i][k].get_data(), u2[k][j].get_data()); // Multiply elements.
                mpz_add(sum, sum, product); // Add to the sum.
                mpz_clear(product); // Clear the product mpz_t to avoid memory leak.
            }
            result[i][j] = fplll::Z_NR<mpz_t>(sum); // Assign the computed sum to the result matrix.
            mpz_clear(sum); // Clear the sum mpz_t to avoid memory leak.
        }
    }

    return result;
}



bool are_matrices_equal(const fplll::ZZ_mat<mpz_t>& matrix1, const fplll::ZZ_mat<mpz_t>& matrix2)
{
    // First, check if the dimensions are the same
    if (matrix1.get_rows() != matrix2.get_rows() || matrix1.get_cols() != matrix2.get_cols()) {
        //return false;
        std::cout<<"size does not match"<<std::endl;
    }

    int rows = matrix1.get_rows();
    int cols = matrix1.get_cols();

    // Iterate through each element and compare
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (mpz_cmp(matrix1[i][j].get_data(), matrix2[i][j].get_data()) != 0) {
                // If any element is not equal, return false
                //std::cout<<i<<","<<j<<std::endl;
                return false;

            }
        }
    }

    // If all elements are equal, return true
    return true;
    //std::cout<<"yeah"<<std::endl;
}

Eigen::MatrixXd zzmat2matrixxd(const fplll::ZZ_mat<mpz_t> &B)
{
    const auto cols = B.get_cols();
    const auto rows = B.get_rows();

    Eigen::MatrixXd R(rows, cols);
    for (auto i = 0; i < rows; i++) // Iterate over rows
    {
        for (auto j = 0; j < cols; j++) // Iterate over columns
        {
            R(i, j) = mpz_get_d((B[i][j]).get_data());
        }
    }

    return R;
}


fplll::ZZ_mat<mpz_t> matrixxd2zzmat(const Eigen::MatrixXd &B)
{
    const auto cols = B.cols();
    const auto rows = B.rows();

    fplll::ZZ_mat<mpz_t> R(rows, cols);
    for (auto i = 0; i < rows; i++) // Corrected loop to iterate over rows
    {
        for (auto j = 0; j < cols; j++) // Corrected loop to iterate over columns
        {
            mpz_t v;
            mpz_init(v);
            // Convert double to long int before using mpz_set_si
            // This is a simplistic conversion; consider handling large or fractional values appropriately
            mpz_set_si(v, static_cast<long>(B(i, j)));
            //mpz_set_si(v, B(i, j));
            R[i][j] = fplll::Z_NR<mpz_t>(v);
            mpz_clear(v); // Clear mpz_t variable to avoid memory leak
        }
    }

    return R;
}

Eigen::MatrixXd lllmatrix2matrixxd(const LLLMatrix &B)
{
    const auto cols = B.cols();
    const auto rows = B.rows();

    Eigen::MatrixXd R(rows, cols);
    for (auto i = 0; i < rows; i++) // Iterate over rows
    {
        for (auto j = 0; j < cols; j++) // Iterate over columns
        {
            R(i, j) = static_cast<double>(B(i, j));
        }
    }

    return R;
}

