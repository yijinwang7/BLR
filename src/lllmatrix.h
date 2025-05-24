// Copyright (c) 2025 Yijin Wang
// This software is released under the MIT License.
// https://opensource.org/licenses/MIT

#ifndef __LLL_MATRIX_H__
#define __LLL_MATRIX_H__

#include <Eigen/Dense>
#include <Eigen/Core>
#include <gmp.h>
#include <memory>

struct ScalarMPZ
{
    mpz_t value;

    ScalarMPZ()
    {
        mpz_init(value);
    }

    ScalarMPZ(const char *str)
    {
        mpz_init_set_str(value, str, 10);
    }

    ScalarMPZ(signed long v)
    {
        mpz_init_set_si(value, v);
    }

    ScalarMPZ(int v)
    {
        mpz_init_set_si(value, v);
    }

    ~ScalarMPZ()
    {
        // mpz_clear(value);
    }

    ScalarMPZ &operator=(const mpz_t &other)
    {
        if (&this->value[0] != &other[0])
        {
            mpz_set(value, other);
        }
        return *this;
    }

    ScalarMPZ &operator=(const ScalarMPZ &other)
    {
        if (this->value != other.value)
        {
            mpz_set(value, other.value);
        }
        return *this;
    }

    ScalarMPZ operator*(const ScalarMPZ &other)
    {
        ScalarMPZ result;

        mpz_mul(result.value, this->value, other.value);
        return result;
    }

    ScalarMPZ operator/(const ScalarMPZ &other)
    {
        ScalarMPZ result;

        mpz_tdiv_q(result.value, this->value, other.value);
        return result;
    }

    ScalarMPZ operator+(const ScalarMPZ &other)
    {
        ScalarMPZ result;

        mpz_add(result.value, this->value, other.value);
        return result;
    }

    ScalarMPZ operator-(const ScalarMPZ &other)
    {
        ScalarMPZ result;

        mpz_sub(result.value, this->value, other.value);
        return result;
    }

    void operator+=(const ScalarMPZ &other)
    {
        ScalarMPZ result;

        mpz_add(result.value, this->value, other.value);
        *this = result;
    }

    void operator-=(const ScalarMPZ &other)
    {
        ScalarMPZ result;

        mpz_sub(result.value, this->value, other.value);
        *this = result;
    }

    bool operator==(const ScalarMPZ &other)
    {
        return mpz_cmp(this->value, other.value) == 0;
    }

    operator double() const
    {
        return mpz_get_d(this->value);
    }

    void get_mpz(mpz_t &v) const
    {
        mpz_init_set(v, this->value);
    }

    friend std::ostream &operator<<(std::ostream &os, const ScalarMPZ &obj)
    {
        size_t len = mpz_sizeinbase(obj.value, 10) + 2;
        auto buf = std::unique_ptr<char[]>(new char[len]);
        mpz_get_str(buf.get(), 10, obj.value);
        os << buf.get();
        return os;
    }


    size_t digits10() const
    {
        size_t minDigits = mpz_sizeinbase(this->value, 10);
        return std::max(minDigits, static_cast<size_t>(mpz_get_d(this->value)));
    }
};

namespace Eigen
{
    template <>
    struct NumTraits<ScalarMPZ>
    {
        typedef ScalarMPZ Real;
        typedef ScalarMPZ NonInteger;
        typedef ScalarMPZ Nested;
        typedef ScalarMPZ Literal;

        enum
        {
            IsInteger = 0,
            IsSigned = 1,
            IsComplex = 0,
            RequireInitialization = 1,
            ReadCost = 1,
            AddCost = 2,
            MulCost = 2,
        };

        static inline int digits10()
        {
            return Real().digits10();
        }
    };
}

inline ScalarMPZ operator*(const ScalarMPZ &lhs, const ScalarMPZ &rhs)
{
    ScalarMPZ result;

    mpz_mul(result.value, lhs.value, rhs.value);
    return result;
}

inline ScalarMPZ operator+(const ScalarMPZ &lhs, const ScalarMPZ &rhs)
{
    ScalarMPZ result;

    mpz_add(result.value, lhs.value, rhs.value);
    return result;
}

inline ScalarMPZ operator-(const ScalarMPZ &lhs, const ScalarMPZ &rhs)
{
    ScalarMPZ result;

    mpz_sub(result.value, lhs.value, rhs.value);
    return result;
}

typedef Eigen::Matrix<ScalarMPZ, Eigen::Dynamic, Eigen::Dynamic> LLLMatrix;
typedef Eigen::Matrix<ScalarMPZ, Eigen::Dynamic, 1> LLLVector;

inline std::ostream& operator<<(std::ostream &os, const LLLMatrix &M) {
    //os << M.rows() << " " << M.cols() << "\n";
    for (int i = 0; i < M.rows(); ++i) {
        for (int j = 0; j < M.cols(); ++j) {
            os << M(i, j)
               << (j + 1 < M.cols() ? ' ' : '\n');
        }
    }
    return os;
}

#endif // __LLL_MATRIX_H__
