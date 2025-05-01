#include <getopt.h>
#include <iostream>
#include <fstream>
#include <string>

#include "BLR.h"
#include "lllmatrix.h"
#include "utils.h"

using namespace std;
using namespace BLR;

static void print_help(const char* prog) {
    cout << "\nUsage: " << prog << " [options] [input_file]\n\n"
         << "Top‐level algorithm (-a, --algorithm):\n"
         << "    lll      L^2\n"
         << "    hlll     householder LLL\n"
         << "    bkz      BKZ reduction\n\n"
         << "When using LLL/HLLL, you can pick an LLL variant:\n"
         << "  -L, --lll-method METHOD   wrapper | heuristic | fast | proved | hkz(not LLL variant, but also a reduction without specifying the b_size)\n\n"
         << "When using BKZ, you can pick a BKZ variant:\n"
         << "  -K, --bkz-method METHOD   default | autoabort | slide | sd\n\n"
         << "General options:\n"
         << "  -d, --delta VAL       Lovász δ (default: 0.99)\n"
         << u8"  -e, --eta   VAL       size‐reduction \u03B7(default: 0.51)\n"
         << "  -t, --theta VAL       θ (for HLLL, default: 0.001)\n"
         << "  -f, --float TYPE      double | dpe | longdouble | mpfr (default: mpfr)\n"
         << "  -p, --precision N     MPFR bits (default: 0)\n"
         << "  -b, --block   N       number of blocks (default: 5)\n"
         << "  -c, --customfre N     custom swap frequency (default: 1 for wrapper variant, 5 otherwise)\n"
         << "  -s, --b_size N        block size, used by BKZ (default: 0)\n"
         << "  -i, --stopCriteria N  stopping‐criterion(default: 0.99 for wrapper variant, 0.999 otherwise)\n"
         << "  -v, --verbose         verbose output\n"
         << "  -h, --help            this message\n\n";
}

LLLMatrix load_matrix(std::istream &in) {
    // 1) peek first non-whitespace char
    char c;
    while (in.get(c) && std::isspace(c)) { }
    if (!in) throw std::runtime_error("Empty input stream");
    in.unget();

    // 2) If it starts with '[', parse bracket-style rows
    if (c == '[') {
        std::vector<std::vector<ScalarMPZ>> rows;
        std::string line;
        while (std::getline(in, line)) {
            // skip empty or non-numeric lines
            if (line.find_first_of("0123456789-+") == std::string::npos)
                continue;

            // trim off leading '[' and trailing ']' or commas
            size_t a = line.find_first_of("0123456789-+");
            size_t b = line.find_last_of("0123456789");
            std::string body = line.substr(a, b - a + 1);

            // split on whitespace
            std::istringstream iss(body);
            std::string tok;
            std::vector<ScalarMPZ> row;
            while (iss >> tok) {
                // use your ScalarMPZ(const char*) constructor
                row.emplace_back(tok.c_str());
            }

            if (!row.empty()) rows.push_back(std::move(row));
        }

        if (rows.empty())
            throw std::runtime_error("No data found in bracketed matrix");

        // all rows must have same length
        size_t m = rows.size();
        size_t n = rows[0].size();
        for (auto &r : rows) {
            if (r.size() != n)
                throw std::runtime_error("Inconsistent row length in matrix");
        }

        // build the Eigen matrix
        LLLMatrix B(m, n);
        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j < n; ++j)
                B(i, j) = rows[i][j];
        return B;
    }

    // 3) Otherwise fall back to “header” style: m n then m×n entries
    size_t m, n;
    if (!(in >> m >> n))
        throw std::runtime_error("Failed to read matrix dimensions");
    LLLMatrix B(m, n);
    std::string tok;
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (!(in >> tok))
                throw std::runtime_error("Unexpected EOF reading matrix entries");
            B(i, j) = ScalarMPZ(tok.c_str());
        }
    }
    return B;
}

//LLLMatrix load_matrix(std::istream &in) {
//    size_t m, n;
//    if (!(in >> m >> n)) {
//        throw std::runtime_error("Failed to read matrix dimensions");
//    }
//
//    LLLMatrix B(m, n);
//    for (size_t i = 0; i < m; ++i) {
//        for (size_t j = 0; j < n; ++j) {
//            std::string token;
//            in >> token;
//            if (!in) throw std::runtime_error("Unexpected EOF reading matrix entries");
//
//            // Option A: convert via raw GMP
//            mpz_t x;
//            mpz_init(x);
//            if (mpz_set_str(x, token.c_str(), 10) < 0) {
//                throw std::runtime_error("Bad integer in input: " + token);
//            }
//            B(i, j) = x;     // uses your ScalarMPZ::operator=(mpz_t&)
//
//            // mpz_clear(x);   // optional, ScalarMPZ assignment makes its own copy
//        }
//    }
//    return B;
//}

std::string to_string(ReductionMethod m) {
    switch (m) {
        case ReductionMethod::LLL:  return "LLL";
        case ReductionMethod::HLLL: return "HLLL";
        case ReductionMethod::BKZ:  return "BKZ";
    }
    return "Unknown";
}

void read_options(int argc, char** argv, Config &cfg, string &input_file) {
    static struct option long_opts[] = {
            {"algorithm",   required_argument, 0, 'a'},
            {"lll-method",  required_argument, 0, 'L'},
            {"bkz-method",  required_argument, 0, 'K'},
            {"delta",       required_argument, 0, 'd'},
            {"eta",         required_argument, 0, 'e'},
            {"theta",       required_argument, 0, 't'},
            {"float",       required_argument, 0, 'f'},
            {"precision",   required_argument, 0, 'p'},
            {"block",       required_argument, 0, 'b'},
            {"b_size",      required_argument, 0, 's'},
            {"customfre",   required_argument, 0, 'c'},
            {"stopCriteria",required_argument, 0, 'i'},
            {"verbose",     no_argument,       0, 'v'},
            {"help",        no_argument,       0, 'h'},
            {0,0,0,0}
    };

    // defaults set by Config constructor
    int opt, idx;
    while ((opt = getopt_long(argc, argv, "a:L:K:d:e:t:f:p:b:i:s:c:vh", long_opts, &idx)) != -1) {
        switch (opt) {
            case 'a': {
                string alg = optarg;
                transform(alg.begin(), alg.end(), alg.begin(), ::tolower);
                if (alg == "lll") {
                    cfg.reductionMethod = ReductionMethod::LLL;
                }
                else if (alg == "hlll") {
                    cfg.reductionMethod = ReductionMethod::HLLL;
                }
                else if (alg == "bkz") {
                    cfg.reductionMethod = ReductionMethod::BKZ;
                }
                else {
                    cerr << "Unknown algorithm: " << optarg << "\n";
                    print_help(argv[0]);
                    exit(1);
                }
                break;
            }

            case 'L': {
                string m = optarg; transform(m.begin(), m.end(), m.begin(), ::tolower);
                if      (m == "wrapper")   cfg.lllMethod = LLLMethod::Wrapper;
                else if (m == "heuristic") cfg.lllMethod = LLLMethod::Heuristic;
                else if (m == "fast")      cfg.lllMethod = LLLMethod::Fast;
                else if (m == "proved")    cfg.lllMethod = LLLMethod::Proved;
                else if (m == "hkz")       cfg.lllMethod = LLLMethod::HKZ;
                else {
                    cerr << "Unknown LLL method: " << optarg << "\n";
                    print_help(argv[0]); exit(1);
                }
                break;
            }

            case 'K': {
                string m = optarg; transform(m.begin(), m.end(), m.begin(), ::tolower);
                if      (m == "default")   cfg.bkzMethod = BKZMethod::Default;
                else if (m == "autoabort") cfg.bkzMethod = BKZMethod::AutoAbort;
                else if (m == "slide")     cfg.bkzMethod = BKZMethod::Slide;
                else if (m == "sd")        cfg.bkzMethod = BKZMethod::SD;
                else {
                    cerr << "Unknown BKZ method: " << optarg << "\n";
                    print_help(argv[0]); exit(1);
                }
                break;
            }

            case 'd': cfg.delta     = stod(optarg); break;
            case 'e': cfg.eta       = stod(optarg); break;
            case 't': cfg.theta     = stod(optarg); break;
            case 'f':
                if      (!strcmp(optarg,"mpfr"))       cfg.ftType = FTType::MPFR;
                else if (!strcmp(optarg,"double"))     cfg.ftType = FTType::Double;
                else if (!strcmp(optarg,"longdouble")) cfg.ftType = FTType::LongDouble;
                else if (!strcmp(optarg,"dpe"))        cfg.ftType = FTType::DPE;
                else if (!strcmp(optarg,"dd"))         cfg.ftType = FTType::DD;
                else if (!strcmp(optarg,"qd"))         cfg.ftType = FTType::QD;
                else {
                    cerr<<"Unknown float type: "<<optarg<<"\n"; print_help(argv[0]); exit(1);
                }
                break;
            case 'p': cfg.precision = stoi(optarg); break;
            case 'b': cfg.numBlocks = stoi(optarg); break;
            case 's': cfg.b_size = stoi(optarg); break;
            case 'c': cfg.customFre = stoi(optarg); break;
            case 'i': cfg.stopCriteria = stoi(optarg); break;
            case 'v': cfg.verbose   = false;         break;
            case 'h': print_help(argv[0]);          exit(0);
            default:  print_help(argv[0]);          exit(1);
        }
    }

    if (optind < argc) input_file = argv[optind];
}

int main(int argc, char **argv) {
    Config config;
    string input_path;
    read_options(argc, argv, config, input_path);

    // Open input
    istream* in = &cin;
    if (!input_path.empty()) {
        static ifstream fin;
        fin.open(input_path);
        if (!fin) {
            cerr << "Error opening file: " << input_path << "\n";
            return 1;
        }
        in = &fin;
    }

    // Load matrix (you must implement this function in BLR)
    LLLMatrix B, raw;
    try {
        //B = load_matrix(*in);
        // load & transpose so that each column is a basis vector
        raw = load_matrix(*in);
        B   = raw.transpose();
    } catch (exception& e) {
        cerr << "Failed to load matrix: " << e.what() << "\n";
        return 1;
    }

    // Run BLR
    if (config.verbose) {
        cerr << "Running with:\n"
             << "  algorithm = " << to_string(config.reductionMethod) << "\n"
             << "  delta     = " << config.delta     << "\n"
             << "  eta       = " << config.eta       << "\n"
             << "  precision = " << config.precision << " bits\n"
             << "  blocks    = " << config.numBlocks << "\n"
             << "  customFre = " << config.customFre << "\n"
             << endl;
    }

    LLLMatrix R;
    try {
        R = blockBLRReduction(B, config);
    } catch (exception& e) {
        cerr << "Reduction error: " << e.what() << "\n";
        return 1;
    }

    // Write result to stdout
    cout << R.transpose() << "\n";


//    //test if sucessfully-reduced
     // Inverse transpose of the transformation matrix
//
//
//    fplll::MatGSO<fplll::Z_NR<mpz_t>, fplll::FP_NR<mpfr_t>> M(fplllB1, arg_u, arg_uinv_t, 0);
//    // one on success
//    int ok = fplll::is_lll_reduced(M, fplll::LLL_DEF_DELTA, fplll::LLL_DEF_ETA);
//
//    std::cout << "LLL-reduced? " << (ok? "yes\n":"no\n");
//
//    // Test if the matrix reduced by fplll is hLLL reduced
//    fplll::ZZ_mat<mpz_t> arg_u1(fplllB1.get_rows(), fplllB1.get_rows()); // Transformation matrix
//    fplll::ZZ_mat<mpz_t> arg_uinv_t1(fplllB1.get_rows(),
//                                     fplllB1.get_rows()); // Inverse transpose of the transformation matrix
//
//
//    // Initialize the MatGSO object with the matrices and flags
//    fplll::MatHouseholder<fplll::Z_NR<mpz_t>, fplll::FP_NR<mpfr_t>> M1(fplllB1, arg_u1, arg_uinv_t1, 0);
//    int isReduced1 = fplll::is_hlll_reduced(M1, fplll::LLL_DEF_DELTA, fplll::LLL_DEF_ETA, fplll::HLLL_DEF_THETA);
//
//    if (isReduced1 != fplll::RED_SUCCESS) {
//        std::cout << "The matrix is not h-LLL reduced by fplll." << std::endl;
//    } else {
//        std::cout << "The matrix is h-LLL reduced by fplll." << std::endl;
//    }
    // 3a) Test LLL‐reducedness (only valid if config.reductionMethod == LLL)
    if (config.reductionMethod == BLR::ReductionMethod::LLL) {
        fplll::ZZ_mat<mpz_t> fplllB1 = lllmatrix2zzmat(R.transpose());
        fplll::ZZ_mat<mpz_t> arg_u(fplllB1.get_rows(), fplllB1.get_rows()); // Transformation matrix
        fplll::ZZ_mat<mpz_t> arg_uinv_t(fplllB1.get_rows(),
                                        fplllB1.get_rows());
        fplll::MatGSO<fplll::Z_NR<mpz_t>, fplll::FP_NR<mpfr_t>>
                M(fplllB1, arg_u, arg_uinv_t, /*verbose=*/0);
        bool ok = fplll::is_lll_reduced(
                M,
                fplll::LLL_DEF_DELTA,    // or config.delta
                fplll::LLL_DEF_ETA       // or config.eta
        );
        std::cout << "LLL-reduced? " << (ok ? "yes\n" : "no\n");
    }

// 3b) Test HLLL‐reducedness (if config.reductionMethod == HLLL)
    if (config.reductionMethod == BLR::ReductionMethod::HLLL) {
        fplll::ZZ_mat<mpz_t> fplllB1 = lllmatrix2zzmat(R.transpose());
        fplll::ZZ_mat<mpz_t> arg_u(fplllB1.get_rows(), fplllB1.get_rows()); // Transformation matrix
        fplll::ZZ_mat<mpz_t> arg_uinv_t(fplllB1.get_rows(),
                                        fplllB1.get_rows());
        fplll::MatHouseholder<fplll::Z_NR<mpz_t>, fplll::FP_NR<mpfr_t>>
                M(fplllB1, arg_u, arg_uinv_t, /*verbose=*/0);
        int hstat = fplll::is_hlll_reduced(
                M,
                fplll::LLL_DEF_DELTA,
                fplll::LLL_DEF_ETA,
                fplll::HLLL_DEF_THETA
        );
        if (hstat == fplll::RED_SUCCESS)
            std::cout << "HLLL-reduced? yes\n";
        else
            std::cout << "HLLL-reduced? no\n";
    }

    return 0;
}



