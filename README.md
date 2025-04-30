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
  - `-t, --theta` θ parameter for HLLL (default: `0.001`)  
  - `-f, --float` choose **`{double, longdouble, mpfr}`** (default: `mpfr`)  
  - `-p, --precision` default: `0`, needs to be specified if used `mpfr`
  - `-b, --block` number of blocks (default: `5`)  
  - `-c, --customfre` custom swap frequency (default: `1` for wrapper variant, `5` otherwise)
  - `-s, --b_size` block size for BKZ (default: `0`)  
  - `-i, --stop_criterion` stopping‐criterion (ratio threshold, default: `0.99` for wrapper variant, `0.999` otherwise)
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

---
## Usage Examples

> **Note:** the **`wrapper`** method (LLL or HLLL) automatically selects an appropriate precision—**omit** `-f` and `-p` when using `wrapper`.

---

# 1. Fast LLL in double-precision

BLR -a lll \
    -L fast \
    -f double \
    -p 0 \
    matrix.txt


# 2. Proved HLLL in MPFR (256-bit)

BLR -a hlll \
    -L proved \
    -f mpfr \
    -p 256 \
    matrix.txt

# 3. wrapper LLL in MPFR (256-bit)
BLR -a lll \
    -L wrapper \
    matrix.txt

# 4. BKZ Reduction (autoabort, block-size = 10)
BLR -a bkz \
    -K autoabort \
    -s 10 \
    matrix.txt

# 5. Piping from latticegen
latticegen u 40 30 \
  | BLR -a lll -L wrapper \
  > reduced_basis.txt


