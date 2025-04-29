# BLRLibrary

A **block lattice-reduction** library and command-line tool built on top of [fplll](https://github.com/fplll/fplll).  
Supports LLL, HLLL, and BKZ variants with customizable floating-point backends and block sizes.

---

## Features

- **CLI front-end** (`BLR`):  
  - `-a, --algorithm` choose **`{lll, hlll, bkz}`**  
  - `-L, --lll-method` choose **`{wrapper, heuristic, fast, proved, hkz}`**  
  - `-K, --bkz-method` choose **`{default, autoabort, slide, sd}`**  
  - `-d, --delta` δ parameter (default: `0.99`)  
  - `-e, --eta` η parameter (default: `0.51`)  
  - `-t, --theta` θ parameter for BKZ/HKZ (default: `0.001`)  
  - `-f, --float` choose **`{double, longdouble, mpfr}`** (default: `mpfr`)  
  - `-p, --precision` MPFR bits (default: `64`)  
  - `-b, --block` number of blocks (default: `5`)  
  - `-c, --customfre` custom swap frequency (default: `0`)  
  - `-v, --verbose`  
  - `-h, --help`

- **Library API** (`include/BLR.h`):  
  - `BLR::Config` — configure your reduction  
  - `BLR::blockBLRReduction(...)` — reduce any `LLLMatrix`

- **Flexible I/O**:  
  - Header format (`m n` + raw entries)  
  - fplll bracketed dumps (`[a b c …]` per line)

---

## Prerequisites

- **C++11** (or later) compiler  
- [CMake](https://cmake.org/) ≥ 3.1  
- [Eigen](https://eigen.tuxfamily.org/)  
- [GMP](https://gmplib.org/) & [MPFR](https://www.mpfr.org/)  
- [fplll](https://github.com/fplll/fplll) (headers & library)  
- OpenMP (optional, for parallel routines)

---

## Building & Installation

```bash
# 1. Clone the repo
git clone https://github.com/yijinwang7/BLRLibrary.git
cd BLRLibrary

# 2. Create & enter a build folder
mkdir build && cd build

# 3. Configure & compile
cmake ..
make 

# 4. (Optional) Install to /usr/local
sudo make install
