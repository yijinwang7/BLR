#include "blockLLLHy_standard.h"

int performReduction(fplll::ZZ_mat<mpz_t> &blocks, const BLR::Config &config, bool useHlll) {
    int status = 0;
    if (!useHlll) {
        // Standard LLL reduction.
        if (config.lllMethod == BLR::LLLMethod::Fast) {
            status = lll_reduction(blocks, config.delta, config.eta,
                                   fplll::LM_FAST, config.getFplllFT(), config.precision, fplll::LLL_DEFAULT);
        } else if (config.lllMethod == BLR::LLLMethod::Wrapper) {
            status = lll_reduction(blocks, config.delta, config.eta,
                                   fplll::LM_WRAPPER, fplll::FT_DEFAULT, 0, fplll::LLL_DEFAULT);
        } else if (config.lllMethod == BLR::LLLMethod::Proved) {
            status = lll_reduction(blocks, config.delta, config.eta,
                                   fplll::LM_PROVED, config.getFplllFT(), config.precision, fplll::LLL_DEFAULT);
        } else if (config.lllMethod == BLR::LLLMethod::Heuristic) {
            status = lll_reduction(blocks, config.delta, config.eta,
                                   fplll::LM_HEURISTIC, config.getFplllFT(), config.precision, fplll::LLL_DEFAULT);
        } else if (config.lllMethod == BLR::LLLMethod::HKZ) {
            // HKZ reduction doesn't require the extra parameters (like b_size)
            status = hkz_reduction(blocks, fplll::HKZ_DEFAULT, config.getFplllFT(), config.precision);
        }
    } else {
//        std::cerr << "[debug] performReduction: HLLL path, method="
//                  << static_cast<int>(config.lllMethod)
//                  << ", delta=" << config.delta
//                  << ", eta="   << config.eta
//                  << ", theta=" << config.theta
//                  << "\n";
        // HLLL reduction.
        if (config.lllMethod == BLR::LLLMethod::Fast) {
            status = hlll_reduction(blocks, config.delta, config.eta,
                                    config.theta, fplll::HLLL_DEF_C,
                                    fplll::LM_FAST, config.getFplllFT(), config.precision, fplll::LLL_DEFAULT, false);
        } else if (config.lllMethod == BLR::LLLMethod::Wrapper) {
            status = hlll_reduction(blocks, config.delta, config.eta,
                                    config.theta, fplll::HLLL_DEF_C,
                                    fplll::LM_WRAPPER, fplll::FT_DEFAULT, 0, fplll::LLL_DEFAULT, false);
        } else if (config.lllMethod == BLR::LLLMethod::Proved) {
            status = hlll_reduction(blocks, config.delta, config.eta,
                                    config.theta, fplll::HLLL_DEF_C,
                                    fplll::LM_PROVED, config.getFplllFT(), config.precision, fplll::LLL_DEFAULT, false);
        }
    }
    return status;
}

LLLMatrix BLR::blockLLLHy_standard(LLLMatrix B, int p, int fre, double stopC, const BLR::Config &config)
{

    int num_block = 2 * p;
    int m = B.cols();
    int n = B.rows();
    int block_size = m / num_block;
    //int sub_size = b_size / num_block;

    std::vector<int> block_sizes(num_block, block_size);  // Initial block size
    int remainder = m % num_block;  // Extra columns

    // Distribute remainder columns
    for (int i = 0; i < remainder; ++i) {
        block_sizes[i]++;
    }

    // Compute starting indices for each block
    std::vector<int> block_starts(num_block, 0);
    for (int i = 1; i < num_block; ++i) {
        block_starts[i] = block_starts[i - 1] + block_sizes[i - 1];
    }

    std::random_device rd;  // Generate a random seed once
    std::mt19937 g(rd());   // Use Mersenne Twister as random engine

    std::vector<int> tags(num_block);
    for (int i = 0; i < num_block; ++i)
    {
        tags[i] = i;
    }

    std::vector<int> new_tags = tags;

    int temp = tags[num_block-1];

    // Shift the second element of each block to the next block
    for (int i = 1; i < num_block-2 ; i = i+2) {
        new_tags[i+2] = tags[i];
    }

    // Place the stored element from the last block into the first block
    new_tags[1] = temp;


    int iter = 0;
    Eigen::MatrixXd B_double = lllmatrix2matrixxd(B);
    Eigen::VectorXd norms = B_double.colwise().norm();
    double od_org = norms.array().log().sum();
    double od;
    std::vector<double> cell_norms(p, 0);
    LLLMatrix B_org(n,m);
    LLLMatrix concatenated, block, block2;
    int start1, size1, start2, size2, current_col, status;


    while (true)
    {
        B_org = B;
        current_col = 0;
        for (int i = 0; i < p; ++i)
        {


            start1 = block_starts[2 * i];
            size1 = block_sizes[2 * i];

            start2 = block_starts[2 * i + 1];
            size2 = block_sizes[2 * i + 1];

            block = B_org.block(0, start1, n, size1);  // First block
            block2 = B_org.block(0, start2, n, size2);  // Second block

            concatenated.resize(n, size1 + size2);
            concatenated << block, block2;


            auto blocks = lllmatrix2zzmat(concatenated.transpose());

            //here, run multiple func
            //status = lll_reduction(blocks, fplll::LLL_DEF_DELTA, fplll::LLL_DEF_ETA, fplll::LM_FAST, fplll::FT_DOUBLE, 0, fplll::LLL_DEFAULT);
            if (config.reductionMethod == BLR::ReductionMethod::LLL) {
                status = performReduction(blocks, config, false);
            } else if (config.reductionMethod == BLR::ReductionMethod::HLLL) {
                status = performReduction(blocks, config, true);
            }
            if (status != fplll::RED_SUCCESS) {
                std::cout << "Reduction failed during block preprocessing: "
                          << fplll::get_red_status_str(status) << std::endl;
            }

            //bkz_reduction(blocks, sub_size, fplll::BKZ_DEFAULT | fplll::BKZ_AUTO_ABORT, fplll::FT_DEFAULT, 0);
            //hkz_reduction(blocks, fplll::HKZ_DEFAULT, fplll::FT_DOUBLE, fplll::FT_DEFAULT);
            //hlll_reduction(blocks, fplll::LLL_DEF_DELTA, fplll::LLL_DEF_ETA, fplll::HLLL_DEF_THETA,fplll::HLLL_DEF_C, fplll::LM_PROVED, fplll::FT_MPFR, 0, fplll::LLL_DEFAULT,false);
            //auto block_reduced = zzmat2lllmatrix(blocks).transpose();
            auto block_reduced = zzmat2lllmatrix(blocks);


            auto block_reduced_trans = block_reduced.transpose();
            Eigen::MatrixXd Bkk = lllmatrix2matrixxd(block_reduced_trans);
            Eigen::VectorXd norms2 = Bkk.colwise().norm();
            cell_norms[i] = norms2.array().log().sum();
            // Place concatenated block back
            B.block(0, current_col, B.rows(), size1 + size2) = block_reduced_trans;
            current_col += size1 + size2;

        }


        iter++;
        if(iter == 1) {
            tags = new_tags;

            // Reorder block_sizes based on tags
            //tags1 = CO2_ordering_new(tags1,num_block);
            std::vector<int> new_block_sizes(num_block);
            std::vector<int> new_block_starts(num_block, 0);
            for (int i = 0; i < num_block; ++i) {
                new_block_sizes[i] = block_sizes[tags[i]];
                if (i > 0) {
                    new_block_starts[i] = block_starts[tags[i]];
                }

            }
            block_sizes = new_block_sizes;  // Update block_sizes to reordered version
            block_starts = new_block_starts;
        }



        od = std::accumulate(cell_norms.begin(), cell_norms.end(), 0.0);

        if ((od / od_org) > stopC) {
            break;
        }
        od_org = od;


        if (iter == 40)
        {
            break;
        }

        if(iter % fre == 0){
            randomPermutation2(B,g);
        }



    }


    LLLMatrix B_trans = B.transpose();
    auto fplllB = lllmatrix2zzmat(B_trans);

    //here also, have options to run multiple func
    if (config.reductionMethod == BLR::ReductionMethod::LLL) {
        status = performReduction(fplllB, config, false);
    } else if (config.reductionMethod == BLR::ReductionMethod::HLLL) {
        status = performReduction(fplllB, config, true);
    }
    if (status != fplll::RED_SUCCESS) {
        std::cout << "Reduction failed during block preprocessing: "
                  << fplll::get_red_status_str(status) << std::endl;
    }
    //status = lll_reduction(fplllB, fplll::LLL_DEF_DELTA, fplll::LLL_DEF_ETA, fplll::LM_FAST, fplll::FT_DOUBLE, 0, fplll::LLL_DEFAULT);
    //bkz_reduction(fplllB, b_size, fplll::BKZ_DEFAULT | fplll::BKZ_AUTO_ABORT, fplll::FT_DEFAULT, 0);
    //hkz_reduction(fplllB, fplll::HKZ_DEFAULT, fplll::FT_DOUBLE, fplll::FT_DEFAULT);
    //hlll_reduction(fplllB, fplll::LLL_DEF_DELTA, fplll::LLL_DEF_ETA, fplll::HLLL_DEF_THETA,fplll::HLLL_DEF_C, fplll::LM_PROVED, fplll::FT_MPFR, 0, fplll::LLL_DEFAULT,false);
    if (status != fplll::RED_SUCCESS)
    {
        std::cout << "LLL reduction failed with error (as the final reduction)" << fplll::get_red_status_str(status) << std::endl;
    }
    LLLMatrix final_B = zzmat2lllmatrix(fplllB).transpose();
    //std::cout <<"iter from block: "<<iter<<std::endl;
    return final_B;

}
