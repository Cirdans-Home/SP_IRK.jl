module SP_IRK
    # This module implements the SP-IRK method for ODEs with implicit
    # fully implicit Runge-Kutta methods in a three step fashion:
    # 1. First a fully parallel diagonalization of a perturbation of the stage
    #    matrix is performed.
    # 2. Then the resulting block-diagonal linear system is solved in parallel
    #    using a direct method.
    # 3. Finally, the solution is used to update the stage values by solving
    #    a low-rank correction problem in the form of a Sylvester matrix equation.
    
    using LinearAlgebra
    using SparseArrays
    using Base.Threads
    using RungeKutta
    using FastGaussQuadrature
    using QuadGK
    using MatrixEquations

    include("schemes.jl")
    include("krylov.jl")
    include("solver.jl")

    export generate_rkcoeff, kpik_sylv_oneside_prealloc, proj_sylvesterc_block_prealloc, proj_sylvesterc_prealloc, rk_lin_thr_solve_mex, rk_nlin_thr_solve_mex
end
