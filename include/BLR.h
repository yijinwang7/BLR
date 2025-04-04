#ifndef BLR_H
#define BLR_H

#include <string>
#include "lllmatrix.h"
#include <fplll.h>

namespace BLR {

    // Define an enum for reduction methods
    enum class ReductionMethod {
        LLL,
        HLLL,
        BKZ
    };

    enum class LLLMethod {
        Wrapper,    // Uses fre = 1
        Heuristic,  // Uses fre = numBlocks
        Fast,       // Uses fre = numBlocks
        Proved,     // Uses fre = numBlocks
        HKZ         //Uses HKZ reduction (doesn't require b_size, so we include HKZ inside the group of LLLMethod)
    };

    enum class BKZMethod {
        Default,    // fplll's default BKZ
        //HKZ,
        AutoAbort,
        Slide,
        SD
    };

    enum class FTType {
        Double,
        MPFR,
        LongDouble,
        DPE,
        DD,
        QD,
        Default
    };

    struct Config {
        // fplll general parameters
        ReductionMethod reductionMethod;
        int b_size;
        //std::string fp_type;
        int precision;
        FTType ftType;

        // BLR-specific parameters
        int numBlocks;
        double stopCriteria;
        //int fre;
        int customFre;      // If user sets this to a nonzero value, use it.

        // Extended options for different reduction variants:
        LLLMethod lllMethod;     // For LLL/HLLL reductions
        BKZMethod bkzMethod;     // For BKZ-related reductions

        // Constructor with defaults:
        Config()
                : reductionMethod(ReductionMethod::LLL),
                  b_size(0), //b_size > 0 implies use of the BKZ-integrated version.
                  //fp_type("double"),
                  ftType(FTType::MPFR),
                  precision(0), //only need to be specified if use FTType MPFR
                  numBlocks(5),
                  stopCriteria(0.999),
                  customFre(0),
                  //fre(1),
                  lllMethod(LLLMethod::Wrapper),
                  bkzMethod(BKZMethod::Default)
        {}
        // Helper method to compute effective fre parameter.
        int getFre() const {
            // If the user provided a custom value, use it.
            if (customFre > 0)
                return customFre;
            // Otherwise, if we're in LLL mode, use 1 for the Wrapper variant.
            if (reductionMethod == ReductionMethod::LLL)
                return (lllMethod == LLLMethod::Wrapper) ? 1 : numBlocks;
                // For BKZ mode, we always use numBlocks.
            else if (reductionMethod == ReductionMethod::BKZ)
                return numBlocks;
            // Default fallback.
            return numBlocks;
        }
        double getStopCriteria() const{
            if (reductionMethod == ReductionMethod::LLL)
                return (lllMethod == LLLMethod::Wrapper) ? 0.99 : stopCriteria;
                // For BKZ mode, we always use 0.999 as default.
            else if (reductionMethod == ReductionMethod::BKZ)
                return stopCriteria;
            else
                return stopCriteria;
        }
        fplll::FloatType getFplllFT() const {
            switch(ftType) {
                case FTType::Double: return fplll::FT_DOUBLE;
                case FTType::MPFR:   return fplll::FT_MPFR;
                case FTType::LongDouble: return fplll::FT_LONG_DOUBLE;
                case FTType::DPE: return fplll::FT_DPE;
                case FTType::DD: return fplll::FT_DD;
                case FTType::QD: return fplll::FT_QD;
                case FTType::Default: return fplll::FT_DEFAULT;
                default: return fplll::FT_MPFR;
            }
        }

    };

    // Declaration of the unified BLR reduction function.
    // The user will pass a lattice matrix and a configuration object.
    LLLMatrix blockBLRReduction(const LLLMatrix &B, const Config &config);

}

#endif // BLR_H
