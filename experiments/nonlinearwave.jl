using SP_IRK
using LinearAlgebra
using SparseArrays
using Base.Threads
# Print animation of the solution
using Plots 
using LaTeXStrings

println("Running Non Linear Wave sample with: ", ARGS, " arguments, nthreads: ", nthreads())

plot_solution = false

method = ARGS[1]
# Convert the ARGS[1] to an integer
s = parse(Int64, ARGS[2])
println("Method: ", method)
println("s: ", s)

dx = 1.0/1023
beta = 10.0

# Discretize 
x = -0.5:dx:0.5
tspan = [0.0, 1.0]
nx = length(x)
nt = 20*nx
println("Number of time steps: ", nt)
println("Number of space steps: ", nx)

# Build 1D Laplacian with Neumann boundary conditions
em = ones(nx)
A  = spdiagm(-1 => em[1:nx-1], 0 => -2.0*em , 1 => em[1:nx-1])
# A[1,1] = -1
# A[nx,nx] = -1
A[1,1] = 1
A[1,2] = 0
A[nx,nx-1] = 0
A[nx,nx] = 1
A = A/dx^2
O = spzeros(nx,nx)
I = spdiagm(0 => ones(nx))
B = [O I; A O]
y0 = [exp.(-100.0.*(x .^ 2)); zeros(nx)]

function odefun_vec(y,p,t,beta,B,nx)
    # Compute the right-hand side of the ODE
    # y = [u; u']
    u = y[1:nx]
    f = B*y + [zeros(nx); beta.*u.^2]
    f[1] = 0.0
    f[nx] = 0.0
    return f
end
odefun(y,p,t) = odefun_vec(y,p,t,beta,B,nx)

# Build the mass matrix
M = spdiagm(0 => ones(2*nx))

# Build the Jacobian function
function jfun(y,p,t,beta,B,nx)
    u = y[1:nx]
    v = [zeros(nx);2*beta.*u]
    J = B + spdiagm(0 => v)
    return J
end

Theta(t,y) = -odefun(y,0,t)
JTheta(t,y) = jfun(y,0,t,beta,B,nx)

t,y,itersyl,ressyl = rk_nlin_thr_solve_mex(method,s,M,Theta,JTheta,tspan,nt,y0)
@btime t,y,itersyl,ressyl = rk_nlin_thr_solve_mex(method,s,M,Theta,JTheta,tspan,nt,y0)

for i in eachindex(itersyl)
    println("Iterations: ", itersyl[i], " Residual: ", ressyl[i])
end

if plot_solution
    println("Plotting solution")
    # Plot the solution
    println("Size of x", size(x))
    println("Size of y", size(y))
    println("Number of time steps: ", nt)
    println("Number of space steps: ", nx)
    anim = @animate for i in 1:nt
        plot(x,y[1:nx,i],ylim=(-0.5,1.0),title = L"u(x,t) \text{ at } t = $(sol.t[i])",xlabel = L"x",ylabel = L"u(x,t)",legend=false)
    end
    gif(anim, "wave.gif", fps = 10)
end
