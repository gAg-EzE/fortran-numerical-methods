# Numerical Methods in Modern Fortran

Three small numerical-analysis programs written in Fortran 2008, each with built-in self-checks:

1. **Bisection**: roots of a polynomial and of a transcendental equation.
2. **Exponential from its Taylor series**, with a recursive factorial, and two ways to compute `cosh(1)`.
3. **Explicit Euler and Heun** for two ODEs with known exact solutions, with a convergence study.

Master 1 Physics (Optoelectronics), Université de Perpignan Via Domitia, Fortran programming course,
by Komlan Gagnon and Kouyarowa Kondjou.

![Convergence of Euler and Heun](figures/ode_convergence.png)

## Build and run

You need `gfortran` and `make`.

```bash
git clone https://github.com/<your-username>/fortran-numerical-methods.git
cd fortran-numerical-methods
make run      # compile (with -Wall -Wextra -std=f2008) and run the three programs
make test     # same as run: each program stops with an error if its checks fail
make plot     # also draw the figures (needs Python with numpy and matplotlib)
make clean
```

## 1. Bisection (`src/bisection.f90`)

`bisection(f, a, b, eps, nmax, n_iter)` takes the function as an argument (through an abstract interface)
and stops when the bracketing interval is shorter than `eps`.

| Equation | Interval | Root found | Iterations |
|---|---|---:|---:|
| x³ − 2.2x = 0 | [0.5, 3] | 1.483243 (exact: √2.2 = 1.483240) | 19 |
| x tan x − 1 = 0 | [2, 4] | 3.425617 | 19 |

## 2. Exponential and factorial (`src/exponential.f90`)

`fact(n)` is a recursive factorial. `expo(x, eps)` sums the Taylor series of exp until the last term is below
`eps` times the sum. Error `expo(x) - exp(x)` for several tolerances:

| x | eps = 1e-1 | eps = 1e-6 | eps = 1e-12 |
|---:|---:|---:|---:|
| −3 | −9.0e-04 | −1.3e-09 | 7.1e-16 |
| −1 | −1.2e-03 | 2.3e-08 | 2.8e-15 |
| 1 | −5.2e-02 | −2.7e-08 | −5.0e-14 |
| 3 | −6.7e-01 | −2.5e-06 | −5.2e-13 |

`cosh(1)` two ways, against the intrinsic function (1.543080634815244):

- `(expo(1) + expo(-1)) / 2` gives 1.543080634815220, which is limited by the series tolerance (error 2e-14);
- the series of cosh, Σ x²ⁱ/(2i)! up to i = 10, gives 1.543080634815244, which matches to all printed digits.

Two numerical points worth knowing:

- **Overflow.** `n!` exceeds the range of a 32-bit integer from n = 13 and of a 64-bit integer from n = 21.
  `fact` uses 64-bit integers and stops with an error beyond 20. `expo` avoids factorials altogether:
  each term is computed from the previous one (`term = term * x / i`).
- **Negative x.** The series of exp(x) for x < 0 alternates in sign, so the terms partly cancel.
  The results above are still accurate for |x| ≤ 3, but for large negative x it is better to compute 1/exp(|x|).

## 3. Euler and Heun (`src/ode.f90`)

Both methods solve u′ = f(t, u) on [0, 1]:

- Euler: uₙ₊₁ = uₙ + h f(tₙ, uₙ)
- Heun: uₙ₊₁ = uₙ + h/2 [f(tₙ, uₙ) + f(tₙ₊₁, uₙ + h f(tₙ, uₙ))]

Test problems:

| Problem | Equation | Exact solution |
|---|---|---|
| 1 | u′ = 1 + 2tu, u(0) = 0 | (√π/2) eᵗ² erf(t) |
| 2 | u′ = −2tu², u(0) = 1 | 1 / (1 + t²) |

The program writes the solution for h = 0.1 and the errors for h = 2⁻ᵏ (k = 1…10) to `results/*.dat`.
The observed order of convergence, computed from the maximum error, tends to the theoretical values:

| Method | Problem 1 | Problem 2 | Theory |
|---|---:|---:|---:|
| Euler | 1.00 | 1.00 | 1 |
| Heun | 2.00 | 2.00 | 2 |

![Numerical and exact solutions for h = 0.1](figures/ode_solutions.png)

## Layout

```text
src/bisection.f90     bisection method and its test
src/exponential.f90   recursive factorial, exp series, cosh
src/ode.f90           Euler and Heun, convergence study
plot_results.py       figures from results/*.dat
results/              output of the programs
figures/              plots
Makefile
```

## License

MIT, see [LICENSE](LICENSE).
