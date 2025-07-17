# GCV - Generalized Cross Validation for Inverse Problems

A MATLAB implementation of Generalized Cross Validation (GCV) methods for parameter selection in 
regularized inverse problems, particularly focused on CT scans restoring and image deblurring applications.

This code is supplementary material for the paper

[Jahn, Kirilin, Convergence of generalized cross-validation with applications to ill-posed integral equations, 2025](https://www.arxiv.org/abs/2506.14558)

A detailed description of all the mathematical grounding can be found there.

## Overview

This repository contains MATLAB code for implementing GCV-based parameter selection methods for solving 
ill-posed inverse problems. The primary focus is CT scan signals restoring and image deblurring using Gaussian blur kernels, with 
efficient implementations for handling large-scale problems through singular value decomposition (SVD) 
techniques.

## Features

- **Generalized Cross Validation (GCV)** implementation for automatic parameter selection
- **Efficient SVD-based** computations for large-scale problems
- **Comparative analysis** between GCV and optimal parameter selection
- **Comprehensive testing framework** with multiple noise levels and problem sizes
- **Image deblurring** with Gaussian blur kernels are deeper optimezed by chunked
- processing for handling large matrices

## Main Functions

## Dependencies

This code requires the following MATLAB toolboxes and external libraries:
- **AIRToolsII** - [Algebraic Iterative Reconstruction Methods](https://github.com/jakobsj/AIRToolsII)
- **IRtools** - [Image Restoration Tools](https://github.com/jnagy1/IRtools)

## Installation and Setup

1. Clone the repository:
```bash
git clone https://github.com/mkirilin/GCV.git
cd GCV
```

2. Set up dependencies in MATLAB:
```matlab
addpath('/path/to/AIRToolsII');
AIRToolsII_setup();
addpath('/path/to/IRtools');
IRtools_setup();
```

3. Run the test suite:
```matlab
run('Tests.m');
```

## Usage Example

```matlab
% Basic GCV parameter selection
[Xgcv, Xopt, err_gcv, err_opt, k_gcv, k_opt] = gcv(S, V, x, bn, m0, m, allSV, coeffs_all);

% The function returns:
% Xgcv    - Solution using GCV-selected parameter
% Xopt    - Solution using optimal parameter
% err_gcv - Relative error for GCV solution
% err_opt - Relative error for optimal solution
% k_gcv   - GCV-selected parameter value
% k_opt   - Optimal parameter value
```

## High-Performance Computing

The repository includes a SLURM batch script (`gcv_test.sh`) for running computations on HPC clusters:

```bash
sbatch gcv_test.sh
```

## Testing and Benchmarking

## Algorithm Details

The implementation uses:
- **Singular Value Decomposition (SVD)** for efficient computation
- **Chunked processing** to handle memory constraints
- **Vectorized operations** for performance optimization
- **Progressive error computation** for all parameter values
- **Discrete Cosine Transform** matrix use for deblurring problem for speeding up and memory efficience

## Applications

This code is particularly useful for:
- Image deblurring and CT sinogram restoration examples
- Signal denoising
- Solving ill-posed inverse problems with GCV for spectral cut-off estimator
- Parameter selection in regularization methods
- Research in computational inverse problems

## Contributing

Contributions are welcome! Please feel free to submit issues, feature requests, or pull requests.

## License

This project is currently without a specified license. Please contact the repository owner for usage permissions.

## Contact

- **Author**: Tim Jahn
- **Email**: jahn@math.tu-berlin.de

- **Author**: Mikhail Kirilin
- **Repository**: [mkirilin/GCV](https://github.com/mkirilin/GCV)
- **Email**: kirilin@math.tu-berlin.de

---

*Note: This repository is part of ongoing research in computational inverse problems and regularization methods.*
