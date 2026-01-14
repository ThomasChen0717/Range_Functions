# Range Functions Analysis Tool

This project provides a suite of tools for analyzing, comparing, and evaluating different interval arithmetic methods for evaluating polynomial ranges on 2D boxes. It is written in Julia.

## Features

- **Interval Evaluation**: Compute ranges of polynomials over 2D boxes using various methods.
- **Method Comparison**: Compare the efficacy (width of computed intervals) of different methods.
- **Analysis**: Visualize the logarithmic distance (logD) vs. logarithmic radius (logr) to analyze convergence properties.

## Supported Methods

- `Lagrange3`: Lagrange interpolation (order 3).
- `Taylor4`: Taylor expansion (order 4).
- `Taylor3`: Taylor expansion (order 3).
- `Taylor2`: Taylor expansion (order 2).
- `Hermite4`: Hermite interpolation (order 4).

## Prerequisites

This project requires **Julia**. The following Julia packages are required:

- `BenchmarkTools`
- `Plots`
- `DataStructures`
- `Colors`
- `Dates`
- `NLsolve`
- `LinearAlgebra`

They should be installed automatically when you run the script for the first time. 
If you encounter any issues, you can install them manually in the Julia REPL:
```julia
using Pkg
Pkg.add(["BenchmarkTools", "Plots", "DataStructures", "Colors", "Dates", "NLsolve", "LinearAlgebra"])
```

## Usage

The main entry point is `range_funcs.jl`. You can run it directly with Julia or use the provided `Makefile`.

### General Syntax

```bash
julia range_funcs.jl <COMMAND> [ARGUMENTS...] [INPUT_FILES...]
```

If no command is specified, it defaults to `intervals`.
If no input file is specified, it defaults to `input.txt`.

### Commands

#### 1. `intervals`
Tests interval evaluation at specific radii around the midpoint of the box defined in the input file.

```bash
julia range_funcs.jl intervals [input.txt]
```
or 
```bash
make intervals [INPUT=input.txt]
```

#### 2. `analyze`
Performs a `logD` vs `logr` analysis for two specified methods and generates a plot in the `imgs/` directory.

Method Arguments default values:
- `METHOD1`: `Lagrange3`
- `METHOD2`: `Taylor4`

```bash
julia range_funcs.jl analyze <METHOD1> <METHOD2> [input.txt]
```
or 
```bash
make analyze [METHOD1=Lagrange3] [METHOD2=Taylor4] [INPUT=input.txt]
```

#### 3. `compare`
Generates a heatmap comparing the ratio of interval widths between two methods over the domain.

Method Arguments default values:
- `METHOD1`: `Lagrange3`
- `METHOD2`: `Hermite4`

```bash
julia range_funcs.jl compare <METHOD1> <METHOD2> [input.txt]
```
or 
```bash
make compare [METHOD1=Lagrange3] [METHOD2=Hermite4] [INPUT=input.txt]
```

**Example:**
```bash
julia range_funcs.jl compare Lagrange3 Hermite4 input.txt
```

## Input File Format

Input files (e.g., `input.txt`) should contain one polynomial definition per line in the following format:

```text
POLYNOMIAL_STRING, x = [MIN, MAX], y = [MIN, MAX]
```

- Lines starting with `#` are treated as comments and ignored.
- The polynomial string should be a standard algebraic expression.

**Example:**
```text
x^2 + y^2 - 1, x = [-1.5, 1.5], y = [-1.5, 1.5]
50*x^10 + 50*y^10 - 1, x = [-1.1, 1.3], y = [-1, 1.4]
```

## Directory Structure

- `range_funcs.jl`: Main entry script.
- `tests.jl`: Core testing and analysis logic.
- `methods/`: Implementation of different interval methods (Lagrange, Taylor, Hermite).
- `SLP/`: Straight-Line Program implementation for polynomial evaluation.
- `imgs/`: Directory where output plots are saved.

## Acknowledgements

This project was developed as part of the research on interval arithmetic methods for polynomial evaluation. 

Mentors:
- Dr. Yap Chee
- Dr. Kai Hormann 
- Bingwei Zhang
