# SP_IRK: Stage Parallel Implicit Runge-Kutta Methods

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://cirdans-home.github.io/SP_IRK.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cirdans-home.github.io/SP_IRK.jl/dev/)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://cirdans-home.github.io/SP_IRK.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cirdans-home.github.io/SP_IRK.jl/dev/)
[![Build Status](https://github.com/cirdans-home/SP_IRK.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/cirdans-home/SP_IRK.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Build Status](https://app.travis-ci.com/cirdans-home/SP_IRK.jl.svg?branch=main)](https://app.travis-ci.com/cirdans-home/SP_IRK.jl)
[![Coverage](https://codecov.io/gh/cirdans-home/SP_IRK.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/cirdans-home/SP_IRK.jl)

Implicit Runge--Kutta (IRK) methods are highly effective for solving stiff ordinary differential equations (ODEs) 
but can be computationally expensive for large-scale problems due to the need to solve coupled algebraic 
equations at each step. This code improves IRK efficiency by leveraging parallelism to decouple stage 
computations and reduce communication overhead. 

We implement here two IRK families:
- symmetric methods
- collocation methods

To treat nonlinear problems the code adopts a simplified Newton method.

## Collaborators

- Fabio Durastante [:email:](mailto:fabio.durastante@unipi.it)
- Mariarosa Mazza [:email:](mailto:mariarosa.mazza@uniroma2.it)

## Reference

If you end up using some of this routines or the ideas contained here, please cite:
```
F. Durastante, M. Mazza. Stage-Parallel Implicit Runge--Kutta methods via low-rank matrix equation corrections.
```


> [!TIP]
> The folder `toeplitzruns` contains some bash script and information to run the examples 
> contained in the accompanying paper on a cluster node whose resources are managed via
> SLURM. Specifically, it has the configuration for using the Toeplitz cluster at the 
> University of Pisa. They should be easily adaptable to other machines. Further information
> are given in a `README.md` file in the folder.

> [!WARNING]
> This code is not a complete Julia package, it is a way to validate the algorithmic 
> proposal and make the numerical tests of the paper accessible and open.