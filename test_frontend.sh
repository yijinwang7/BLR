#!/usr/bin/env bash
set -euo pipefail

# Path to BLR executable, may be re-adjusted according to the build path
BLR_CMD="./cmake-build-debug/BLR"

# Loop over the six latticegen types
for type in r u q n t s; do
  case "$type" in
    r) LG_ARGS=(r 20 5) ;;             # knapsack, d=20, bits=5
    u) LG_ARGS=(u 20 5) ;;             # uniform, d=20, bits=5
    q) LG_ARGS=(q 30 10 5 b) ;;        # q-ary, d=30,k=10,bits=5, random q
    n) LG_ARGS=(n 20 5 b) ;;             # ntru-like, d=20,k=5
    t) LG_ARGS=(t 20 0.5) ;;           # triangular, d=20, density=0.5
    s) LG_ARGS=(s 20 10 30) ;;         # simdioph, d=20,r=10,bitlen=30
  esac

  echo
  echo "======================================================"
  echo " Testing lattice type '$type' → latticegen ${LG_ARGS[*]}"
  echo "======================================================"

  # 1) LLL variants
  for method in hkz wrapper heuristic fast proved; do
    echo
    echo "→ LLL / $method"
    if [ "$method" = "wrapper" ]; then
      # wrapper → no -f or -p
      latticegen "${LG_ARGS[@]}" \
        | $BLR_CMD -a lll -L wrapper \
        || echo "!! Error in LLL/wrapper for type $type"
    elif [ "$method" = "hkz" ]; then
      # pure HKZ reduction → pick a float backend, here double + prec=0
      latticegen "${LG_ARGS[@]}" \
        | $BLR_CMD \
            -a lll \
            -L hkz \
            -f double \
            -p 0 \
        || echo "!! Error in LLL/hkz for type $type"
    else
      # pick float and precision for the non-wrapper methods
      case $method in
        heuristic) FLOAT=dpe    ; PREC=0 ;;
        fast)      FLOAT=double ; PREC=0 ;;
        proved)    FLOAT=mpfr   ; PREC=0 ;;
      esac

      latticegen "${LG_ARGS[@]}" \
        | $BLR_CMD \
            -a lll \
            -L $method \
            -f $FLOAT \
            -p $PREC \
        || echo "!! Error in LLL/$method for type $type"
    fi
  done

  # 2) HLLL variants
  for method in wrapper fast proved; do
    echo
    echo "→ HLLL / $method"
    if [ "$method" = "wrapper" ]; then
      latticegen "${LG_ARGS[@]}" \
        | $BLR_CMD -a hlll -L wrapper \
        || echo "!! Error in HLLL/wrapper for type $type"
    else
      case $method in
        fast)   FLOAT=double ; PREC=0 ;;
        proved) FLOAT=mpfr   ; PREC=0 ;;
      esac

      latticegen "${LG_ARGS[@]}" \
        | $BLR_CMD \
            -a hlll \
            -L $method \
            -f $FLOAT \
            -p $PREC \
        || echo "!! Error in HLLL/$method for type $type"
    fi
  done

  # 3) BKZ with b_size=10
  echo
  echo "→ BKZ variants (b_size=10)"
  for bkz_m in default autoabort slide sd; do
    echo
    echo "---- BKZ / $bkz_m ----"
    latticegen "${LG_ARGS[@]}" \
      | $BLR_CMD \
          -a bkz \
          -K $bkz_m \
          -s 10 \
	  -f mpfr \
	  -p 64 \
      || echo "!! Error in BKZ/$bkz_m for type $type"
  done

done

echo
echo "=== quick_test.sh done ==="

