#ifndef __LLL_REDUCTION_UTILS_H__
#define __LLL_REDUCTION_UTILS_H__

#include <Eigen/Dense>
#include <fplll.h>
#include "lllmatrix.h"

LLLMatrix zzmat2lllmatrix(const fplll::ZZ_mat<mpz_t> &B);

fplll::ZZ_mat<mpz_t> lllmatrix2zzmat(const LLLMatrix &B);

fplll::ZZ_mat<mpz_t> multiply_matrices(const fplll::ZZ_mat<mpz_t>& u1, const fplll::ZZ_mat<mpz_t>& u2);

bool are_matrices_equal(const fplll::ZZ_mat<mpz_t>& matrix1, const fplll::ZZ_mat<mpz_t>& matrix2);

Eigen::MatrixXd zzmat2matrixxd(const fplll::ZZ_mat<mpz_t> &B);

fplll::ZZ_mat<mpz_t> matrixxd2zzmat(const Eigen::MatrixXd &B);

Eigen::MatrixXd lllmatrix2matrixxd(const LLLMatrix &B);

#endif // __LLL_REDUCTION_UTILS_H__