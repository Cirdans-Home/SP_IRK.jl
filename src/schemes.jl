"""
Generate the coefficients for an s-stage Runge-Kutta method.
    A, b, c, X, ev = generate_rkcoeff(s,type) # Only for Gauss
    A, b, c, M, B, W = generate_rkcoeff(s,type) # For all the other methods
 
# INPUT:
- s: number of stages
- type: type of Runge-Kutta method, accepted types are: 
    - Gauss (default) "gauss"
    - Radau IA  "radauIA"
    - Radau IIA "radauIIA"
    - Lobatto IIIA "lobattoIIIA"
    - Lobatto IIIB "lobatto IIIB"
    - Lobatto IIIC "lobatto IIIC"
    - Lobatto IIID "lobatto IIID"

# OUTPUT:
- A: entries of the Butcher Tableau
- b: entries of the Butcher Tableau
- c: entries of the Butcher Tableau
- X: eigenvectors of the centrosymmetric matrix associated with the matrix A (Only for Gauss)
- ev: eigenvalues of the centrosymmetric matrix associated with the matrix A (Only for Gauss)
- M: M = W' B A W transformed matrix A of the Butcher tableau
- W: W-transformation matrix
- B: diagonal matrix with the b coefficients of the Butcher tableau
"""
function generate_rkcoeff(s,type="gauss")
    
    D = Diagonal(ones(s))

    # Switch over the type of quadrature rule
    if type == "gauss"
        # Gauss-Legendre collocation Runge-Kutta method
        c, b = gausslegendre(s)
        c = map(x -> (x + 1) / 2, c)
        b = b / 2

        A = zeros(s, s)
        for i = 1:s
            for j = 1:s
                A[i, j], err = quadgk(x -> lagrangePolynomial(c, j, x), 0, c[i], rtol=1e-12, atol=1e-12)
            end
        end
    
        e = ones(s)
        M = A - 0.5 * e * transpose(b)
        ev, X = eigen(M)
        # Set the real part of ev to zeros
        ev = complex.(zeros(s), imag(ev))
    
        return A, b, c, X, ev, D
    elseif type == "radauIA"
        # Gauss-Radau collocation Runge-Kutta method
        tab = TableauRadauIA(s)
        b = tab.b
        c = tab.c
        A = convert(Matrix{Float64},tab.a)
        B,W = wtransform(s,b,c)
        M = convert(Matrix{Float64},W'*B*A*W)
        # Set the entries of M that are of absolute value smaller than 1e-13 to zeros
        M = map(x -> abs(x) < 1e-13 ? 0 : x, M)

        return A, b, c, M, B, W, D

    elseif type == "radauIIA"

        # Gauss-Radau collocation Runge-Kutta method
        tab = TableauRadauIIA(s)
        b = tab.b
        c = tab.c
        A = tab.a
        B,W = wtransform(s,b,c)
        M = W'*B*A*W
        # Set the entries of M that are of absolute value smaller than 1e-13 to zeros
        M = map(x -> abs(x) < 1e-13 ? 0 : x, M)

        return A, b, c, M, B, W, D

    elseif type == "lobattoIIIA"

        # Gauss-Lobatto collocation Runge-Kutta method
        tab = TableauLobattoIIIA(s)
        b = tab.b
        c = tab.c
        A = tab.a
        B,W = wtransform(s,b,c)
        M = W'*B*A*W
        # Set the entries of M that are of absolute value smaller than 1e-13 to zeros
        M = map(x -> abs(x) < 1e-13 ? 0 : x, M)
        D[s,s] = (2.0*s-1)/(s-1.0)

        return A, b, c, M, B, W, D

    elseif type == "lobattoIIIB"

        # Gauss-Lobatto collocation Runge-Kutta method
        tab = TableauLobattoIIIB(s)
        b = tab.b
        c = tab.c
        A = tab.a
        B,W = wtransform(s,b,c)
        M = W'*B*A*W
        # Set the entries of M that are of absolute value smaller than 1e-13 to zeros
        M = map(x -> abs(x) < 1e-13 ? 0 : x, M)
        D[s,s] = (2.0*s-1)/(s-1.0)

        return A, b, c, M, B, W, D

    elseif type == "lobattoIIIC"

        # Gauss-Lobatto collocation Runge-Kutta method
        tab = TableauLobattoIIIC(s)
        b = tab.b
        c = tab.c
        A = tab.a
        B,W = wtransform(s,b,c)
        M = W'*B*A*W
        # Set the entries of M that are of absolute value smaller than 1e-13 to zeros
        M = map(x -> abs(x) < 1e-13 ? 0 : x, M)
        D[s,s] = (2.0*s-1)/(s-1.0)

        return A, b, c, M, B, W, D

    elseif type == "lobattoIIID"

        # Gauss-Lobatto collocation Runge-Kutta method
        tab = TableauLobattoIIID(s)
        b = tab.b
        c = tab.c
        A = tab.a
        B,W = wtransform(s,b,c)
        M = W'*B*A*W
        # Set the entries of M that are of absolute value smaller than 1e-13 to zeros
        M = map(x -> abs(x) < 1e-13 ? 0 : x, M)
        D[s,s] = (2.0*s-1)/(s-1.0)

        return A, b, c, M, B, W, D

    else
        error("Invalid type of quadrature rule")
    end
end

"""
Evaluates the jth Lagrange polynomial at point x.

# Arguments
- c: A vector of nodes
- j: The index of the Lagrange polynomial to evaluate
- x: The point at which to evaluate the Lagrange polynomial
# Outputs:
- Lj: The value of the jth Lagrange polynomial at point x
"""
function lagrangePolynomial(c, j, x) 
    # Number of nodes
    n = length(c)
    
    # Initialize Lj to 1
    Lj = 1
    
    # Loop over all nodes except the jth node
    for m = 1:n
        if m != j
            # Compute the product for the Lagrange basis polynomial
            Lj = Lj * (x - c[m]) / (c[j] - c[m])
        end
    end
    
    return Lj
end

"""
Compute the W-transformation matrix for a given Runge-Kutta method.

#INPUT:
s : Int
    The order of the Runge-Kutta method
b : Array{Float64, 1}
    The b coefficients of the Runge-Kutta method
c : Array{Float64, 1}
    The c coefficients of the Runge-Kutta method

"""
function wtransform(s, b, c)
    # Compute the B matrix as the diagonal matix with elements taken from B
    B = diagm(b)
    # Use the recurrence relation for the normalized Legendre polynomials on the [0,1] interval to compute
    # the matrix P_{j-1}(c[i]) for j = 1,...,s and i = 1,...,s
    W = orthonormal_legendre(c, s)

    return B, W
end

"""
orthonormal_legendre(c, s)

Compute the orthonormal Legendre polynomials.

This function computes the orthonormal Legendre polynomials up to degree `s` using the recurrence relation method.
The Legendre polynomials are a set of orthonormal polynomials defined on the interval [0, 1].
The coefficients `c` are the evaluation points at which the polynomials are computed.

# Arguments
- `c`: An array of evaluation points.
- `s`: The degree of the highest polynomial to compute.

# Returns
- `W`: An array of size `(length(c), s)` containing the computed orthonormal Legendre polynomials.

"""
function orthonormal_legendre(c, s)
    W = zeros(s, s)
    # coefficients of the recurrence relation
    @inline A(n) = (4.0 * n + 2.0) / (n + 1)
    @inline B(n) = -(2.0 * n + 1) / (n + 1)
    @inline C(n) = n / (n + 1)
    # First column is one irrespectively of the evaluation point
    W[:, 1] = ones(s, 1)
    # Second column
    if s >= 2
        W[:, 2] = A(0) * c .+ B(0)
    end
    # Three-terms recurrence
    for j in 2:s-1
        n = j - 1
        W[:, j+1] = (A(n) * c .+ B(n)) .* W[:, j] - C(n) * W[:, j-1]
    end
    # Normalize the columns
    for j in 0:s-1
        W[:, j+1] = W[:, j+1] * sqrt(2 * j + 1)
    end

    return W

end