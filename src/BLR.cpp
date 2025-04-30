#include "BLR.h"
#include "blockLLLHy_BKZ.h"  // Ensure this header declares both versions.
#include "blockLLLHy_standard.h"
#include "randomPermutation.h"
#include "utils.h"

namespace BLR {
    LLLMatrix blockBLRReduction(const LLLMatrix &B, const Config &config) {
        if (config.reductionMethod == ReductionMethod::LLL||config.reductionMethod == ReductionMethod::HLLL) {
            // For LLL/HLLL/HKZ, call the standard version.
            return blockLLLHy_standard(B, config.numBlocks, config.getFre(), config.getStopCriteria(), config);
        } else if (config.reductionMethod == ReductionMethod::BKZ) {
            // For BKZ, call the BKZ version.
            return blockLLLHy_BKZ(B, config.numBlocks, config.getFre(), config.getStopCriteria(), config.b_size, config);
        }
        return B;
    }
}
