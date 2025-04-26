"""
    kpik_sylv_oneside_prealloc(Am,As, B, rhs1, rhs2, m, tol, U, H)

Solve the Sylvester equation of the form `AX + XB = C` using a one-sided Krylov 
projection method with preallocated memory for efficiency.

# Arguments
- `Am::Function`: A function that computes the action of the matrix `A` on a vector.
- `As::Function`: A function that computes the action of the matrix `A^{-1}` on a vector.
- `B::AbstractMatrix`: The matrix `B` in the Sylvester equation.
- `rhs::AbstractVector`: The first low-rank right-hand side of the Sylvester equation.
- `rhs2::AbstractVector`: The second low-rank right-hand side of the Sylvester equation.
- `m::Integer`: The maximum number of Krylov iterations.
- `tol::Float64`: The tolerance for convergence.
- `U::AbstractMatrix`: The Krylov basis matrix, which will be filled during the computation.
- `H::AbstractMatrix`: The matrix that will store the coefficients of the Krylov basis.
# Returns
- `X::AbstractMatrix`: The solution matrix to the Sylvester equation.
- `iter::Integer`: The number of iterations performed.
- `res::Float64`: The residual norm of the solution.

# Notes
This function is designed for performance-critical applications where 
preallocation of memory is necessary to avoid repeated allocations during 
iterative computations. Ensure that the dimensions of `A`, `B`, `C`, and `X` 
are compatible with the Sylvester equation.
"""
function kpik_sylv_oneside_prealloc(Am::Function, As::Function, B::AbstractMatrix, rhs1::AbstractVector, rhs2::AbstractMatrix, m::Integer, tol::Float64, U::AbstractMatrix, H::AbstractMatrix)
    nrmb1 = norm(rhs1)
    nrmb2 = norm(rhs2)
    normb = sqrt(nrmb1 * nrmb2)

    n = size(rhs1, 1)
    sh = size(rhs1, 2)
    sh1 = size(rhs2, 2)

    if sh1 != sh
        error("The two right-hand sides must have the same number of columns")
    end

    odds = []
    res = zeros(m + 2)
    kA_max = m

    s = 2 * sh
    Up = zeros(n, sh + s - 1)

    Is = spdiagm(0 => ones(s))
    rhsA = As(rhs1)
    F = qr(hcat(rhs1, rhsA))
    U[1:n, 1:s] = Matrix(F.Q)
    beta = F.R
    ibeta = beta \ Is
    beta = beta[1:sh, 1:sh]
    L = zeros((m + 1) * s, m * s)
    T = zeros(m * s, m * s)

    beta2 = beta
    rho = 0

    for j = 1:m
        jms = (j - 1) * s + 1
        j1s = (j + 1) * s
        global js = j * s
        js1 = js + 1
        jsh = (j - 1) * s + sh

        # Sequence in A
        Up[1:n, 1:sh] = Am(U[1:n, jms:jsh])
        Up[1:n, sh+1:s] = As(U[1:n, jsh+1:js])
        for l = 1:2
            k_min = max(1, j - kA_max)
            for kk = k_min:j
                k1 = (kk - 1) * s + 1
                k2 = kk * s
                coef = U[1:n, k1:k2]' * Up
                H[k1:k2, jms:js] = H[k1:k2, jms:js] + coef
                Up = Up - U[1:n, k1:k2] * coef
            end
        end
        if j <= m
            F = qr(Up)
            U[1:n, js1:j1s] = Matrix(F.Q)
            H[js1:j1s, jms:js] = Matrix(F.R)
            hinv = inv(H[js1:j1s, jms:js])
        end
        I = spdiagm(0 => ones(js + s))
        if j == 1
            L[1:j*s+sh, (j-1)*sh+1:j*sh] = hcat(H[1:s+sh, 1:sh] / ibeta[1:sh, 1:sh], vcat(spdiagm(0 => ones(sh)), spzeros(s, sh)) / ibeta[1:sh, 1:sh]) * ibeta[1:s, sh+1:s]
        else
            L[1:j*s+s, (j-1)*sh+1:j*sh] = L[1:j*s+s, (j-1)*sh+1:j*sh] + H[1:j*s+s, jms:jms-1+sh] * rho
        end
        odds = Vector(vcat(odds, jms:(jms-1+sh)))
        evens = Vector(1:js)

        deleteat!(evens, odds) # evens[odds] = [] in Matlab

        T[1:js+s, odds] = H[1:js+s, odds] # Odd columns
        T[1:js+sh, evens] = L[1:js+sh, 1:j*sh] # Even columns

        L[1:j*s+s, j*sh+1:(j+1)*sh] = (I[1:j*s+s, (js-sh+1):js] - T[1:js+s, 1:js] * H[1:js, js-sh+1:js]) * hinv[sh+1:s, sh+1:s]

        rho = hinv[1:sh, 1:sh] \ hinv[1:sh, sh+1:s]


        k = j

        Eks = zeros(k * s, sh)
        Eks[1:sh, 1:sh] = diagm(0 => ones(sh))

        global Y = sylvc(T[1:js, 1:js], B, Eks * beta2 * rhs2')

        cc = hcat(H[js1:j1s, js-s+1:js-sh], L[js1:j1s, (j-1)*sh+1:j*sh])

        res[k] = norm(cc * Y[js-s+1:js, :]) / normb

        if res[k] < tol
            global iter = k
            X = U[1:n, 1:js] * Y
            res = res[iter]

            return X, iter, res
        end
    end

    X = U[1:n, 1:js] * Y
    res = res[iter]

    return X, iter, res

end

"""
    proj_sylvesterc_block_prealloc(h, L, M, A, C1, C2, u)

Solve the Sylvester equation:

    h⁻¹ L⁻¹ M X + X Aᵀ = C1 C2ᵀ

with a low-rank right-hand side, where `L` and `M` are large and sparse matrices, 
and `A` is a small and dense matrix.

The solution is computed using a one-sided Krylov projection algorithm with 
`h⁻¹ L⁻¹ M` and `u`, and the Bartels and Stewart algorithm on the projected matrix.

# Arguments
- `h::Float64`: A scalar value used in the equation.
- `M::AbstractMatrix`: A large, sparse matrix.
- `L::AbstractMatrix`: A large, sparse matrix.
- `At::AbstractMatrix`: A small, dense matrix (transpose of `A`).
- `C1::AbstractMatrix`: A matrix representing the first component of the low-rank right-hand side.
- `C2::AbstractMatrix`: A matrix representing the second component of the low-rank right-hand side.
- `maxit::Integer`: The maximum number of Krylov iterations.
- `eps::Float64`: The tolerance for convergence.
- `V::AbstractMatrix`: The Krylov basis matrix, which will be filled during the computation.
- `H::AbstractMatrix`: The matrix that will store the coefficients of the Krylov basis.

# Returns
- The solution matrix `X` to the Sylvester equation.
- The number of iterations performed.
- The residual norm of the solution.

# Notes
This function is optimized for cases where `L` and `M` are large and sparse, 
and the right-hand side has a low-rank structure.
"""
function proj_sylvesterc_block_prealloc(h::Float64, M::AbstractMatrix, L, At::AbstractMatrix, C1::AbstractVector, C2::AbstractMatrix, maxit::Integer, eps::Float64, V::AbstractMatrix, H::AbstractMatrix)
    # Rank of rhs
    rk = size(C1, 2)
    # First Krylov vector is the right-hand side normalized
    normb = norm(C1, 2) * norm(C2, 2)
    QRC1 = qr(C1)
    V[:, 1:rk] = Matrix(QRC1.Q)
    e1 = zeros(rk * maxit, rk)
    e1[1:rk, :] = Matrix(QRC1.R)

    # Loop up to maxit
    for j = 1:maxit-1
        # Apply L to the current Krylov vector
        W = (L \ (M * V[:, ((j-1)*rk+1):j*rk])) / h
        # Orthogonalize the new vector against the previous ones
        for i = 1:j
            H[((i-1)*rk+1):i*rk, ((j-1)*rk+1):j*rk] = transpose(V[:, ((i-1)*rk+1):i*rk]) * W
            W = W - V[:, ((i-1)*rk+1):i*rk] * H[((i-1)*rk+1):i*rk, ((j-1)*rk+1):j*rk]
        end
        # Normalize the new vector
        W = qr(W)
        H[(j*rk+1):(j+1)*rk, ((j-1)*rk+1):j*rk] = Matrix(W.R)
        V[:, (j*rk+1):(j+1)*rk] = Matrix(W.Q)
        # Solve the projected Sylvester equation
        RHS = e1[1:j*rk, :] * transpose(C2)
        Y = sylvc(H[1:j*rk, 1:j*rk], At, RHS)
        # Compute the residual
        res = norm(H[(j*rk+1):(j+1)*rk, ((j-1)*rk+1):j*rk] * Y[((j-1)*rk+1):j*rk, :]) / normb
        # println("Iteration ", j, " Residual ", res)
        # Check for convergence
        if res < eps
            X = V[:, 1:j*rk] * Y
            return X, j, res
        end
    end

    # We did not converge within maxit
    Y = sylvc(H[1:rk*maxit, 1:rk*maxit], transpose(A), e1[1:rk*maxit] * C2')
    X = V[:, 1:rk*maxit] * Y
    res = norm((L \ (M * X)) / h + X * transpose(A) - C1 * C2')

    return X, maxit, res

end
"""
    solve_sylvester(V, H, L, M, A, u, v, h)

Solve the Sylvester equation:

    h⁻¹ L⁻¹ M X + X Aᵀ = uvᵀ

using preallocated memory for `V` and `H`. This function is designed for cases where `L` and `M` 
are large and sparse matrices, `A` is a small and dense matrix, and the right-hand side is rank-1.

The solution is computed using a one-sided Krylov projection algorithm with `h⁻¹ L⁻¹ M` and `u`, 
and the Bartels and Stewart algorithm is applied to the projected matrix.

# Arguments
- `V`: Preallocated matrix to store the Krylov basis.
- `H`: Preallocated matrix to store the Hessenberg matrix from the Krylov projection.
- `L`: Sparse matrix used in the left-hand side of the equation.
- `M`: Sparse matrix used in the left-hand side of the equation.
- `A`: Small dense matrix used in the right-hand side of the equation.
- `u`: Vector defining the rank-1 right-hand side.
- `v`: Vector defining the rank-1 right-hand side.
- `h`: Scalar factor applied to the left-hand side.

# Returns
- `X`: Solution matrix to the Sylvester equation.
- The number of iterations performed.
- The residual norm of the solution.

# Notes
This implementation is optimized for performance by leveraging preallocated memory and efficient 
algorithms for both the Krylov projection and the solution of the reduced Sylvester equation.
"""
function proj_sylvesterc_prealloc(h::Float64, M::AbstractMatrix, L, A::AbstractMatrix, u::AbstractArray, v::AbstractArray, maxit::Integer, eps::Float64, V::AbstractMatrix, H::AbstractMatrix)
    # First Krylov vector is the right-hand side normalized
    normb = sqrt(dot(v, v) * dot(u, u))
    beta = norm(u)
    V[:, 1] = u / beta
    e1 = zeros(maxit, 1)
    e1[1] = beta
    # Loop up to maxit
    for j = 1:maxit-1
        # Apply L to the current Krylov vector
        w = (L \ (M * (V[:, j] / h)))
        # Orthogonalize the new vector against the previous ones
        for i = 1:j
            H[i, j] = dot(V[:, i], w)
            w = w - H[i, j] * V[:, i]
        end
        # Normalize the new vector
        H[j+1, j] = norm(w)
        V[:, j+1] = w / H[j+1, j]
        # Solve the projected Sylvester equation
        Y = sylvc(H[1:j, 1:j], A, e1[1:j] * v')
        # Compute the residual
        res = abs(H[j+1, j]) * norm(Y[j, :], 2) / normb
        # Check for convergence
        if res < eps
            X = V[:, 1:j] * Y
            return X, j, res
        end
    end

    # We did not converge within maxit
    Y = sylvc(H[1:maxit, 1:maxit], A, e1[1:maxit] * v')
    X = V[:, 1:maxit] * Y
    res = norm((L \ (M * X)) / h + X * A - u * v')

    return X, maxit, res

end
"""
    solve_sylvester(V, H, L, M, A, u, v, h)

Solve the Sylvester equation:

    h⁻¹ L⁻¹ M X + X Aᵀ = uvᵀ

using preallocated memory for `V` and `H`. This function is designed for cases where `L` and `M` 
are large and sparse matrices, `A` is a small and dense matrix, and the right-hand side is rank-1.

The solution is computed using a one-sided Krylov projection algorithm with `h⁻¹ L⁻¹ M` and `u`, 
and the Bartels and Stewart algorithm is applied to the projected matrix.

# Arguments
- `V`: Preallocated matrix to store the Krylov basis.
- `H`: Preallocated matrix to store the Hessenberg matrix from the Krylov projection.
- `L`: Sparse matrix used in the left-hand side of the equation.
- `M`: Sparse matrix used in the left-hand side of the equation.
- `A`: Small dense matrix used in the right-hand side of the equation.
- `u`: Vector defining the rank-1 right-hand side.
- `v`: Vector defining the rank-1 right-hand side.
- `h`: Scalar factor applied to the left-hand side.

# Returns
- `X`: Solution matrix to the Sylvester equation.
- The number of iterations performed.
- The residual norm of the solution.

# Notes
This implementation is optimized for performance by leveraging preallocated memory and efficient 
algorithms for both the Krylov projection and the solution of the reduced Sylvester equation.
"""
function proj_sylvesterc(h::Float64, M::AbstractMatrix, L, A::AbstractMatrix, u::AbstractArray, v::AbstractArray, maxit::Integer, eps::Float64)
    # Allocate memory for Krylov basis and Hessenberg matrix
    n = size(M, 1)
    V = zeros(n, maxit)
    H = zeros(maxit+1, maxit)
    # First Krylov vector is the right-hand side normalized
    normb = sqrt(dot(v, v) * dot(u, u))
    beta = norm(u)
    V[:, 1] = u / beta
    e1 = zeros(maxit, 1)
    e1[1] = beta
    # Loop up to maxit
    for j = 1:maxit-1
        # Apply L to the current Krylov vector
        w = (L \ (M * (V[:, j] / h)))
        # Orthogonalize the new vector against the previous ones
        for i = 1:j
            H[i, j] = dot(V[:, i], w)
            w = w - H[i, j] * V[:, i]
        end
        # Normalize the new vector
        H[j+1, j] = norm(w)
        V[:, j+1] = w / H[j+1, j]
        # Solve the projected Sylvester equation
        Y = sylvc(H[1:j, 1:j], A, e1[1:j] * v')
        # Compute the residual
        res = abs(H[j+1, j]) * norm(Y[j, :], 2) / normb
        # Check for convergence
        if res < eps
            X = V[:, 1:j] * Y
            return X, j, res
        end
    end

    # We did not converge within maxit
    Y = sylvc(H[1:maxit, 1:maxit], A, e1[1:maxit] * v')
    X = V[:, 1:maxit] * Y
    res = norm((L \ (M * X)) / h + X * A - u * v')

    return X, maxit, res

end