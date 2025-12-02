# VecCmplx
Vectorized complex multiplication.

(... work in progress ...)

The library should compile on any `AVX512F`-compatible system with either `icx` or `gcc`.

Requires [libpvn](https://github.com/venovako/libpvn) and [SLEEF](https://sleef.org) with quadruple precision support
(please set the ``SLEEF`` and ``sleef=0`` make variables for libpvn).

This work has been supported in part by Croatian Science Foundation under the project IP-2014-09-3670 ([MFBDA](https://web.math.pmf.unizg.hr/mfbda/)).
