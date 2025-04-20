using Gridap
using SP_IRK
using LinearAlgebra
using SparseArrays
using BenchmarkTools
using Base.Threads

println("Running FEM sample with: ",ARGS, " arguments, nthreads: ", nthreads())

method = ARGS[1]
# Convert the ARGS[1] to an integer
s = parse(Int64,ARGS[2])
println("Method: ", method)
println("s: ", s)

# Using Gridap to assemble the problem
model = DiscreteModelFromFile("./experiments/model.json")
order = parse(Int64,ARGS[3])
println("Order: ", order)
reffe = ReferenceFE(lagrangian,Float64,order)
V0 = TestFESpace(model,reffe;conformity=:H1,dirichlet_tags="sides")
g(x) = 2.0
Ug = TrialFESpace(V0,g)
degree = order+1
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)
neumanntags = ["circle", "triangle", "square"]
Γ = BoundaryTriangulation(model,tags=neumanntags)
dΓ = Measure(Γ,degree)
u0f(x) = 2.0
writevtk(model,"model")
global ftrue(x,t) = (1.5-x[1]).*(1-x[2]).*(1-x[3]) + 0.5*sin(2*pi*t)
f(x) = ftrue(x,0.0)
h(x) = 3.0
m(u,v) = ∫( v*u )*dΩ
a(u,v) = ∫( ∇(v)⋅∇(u) )*dΩ
b(v) = ∫( v*f )*dΩ + ∫( v*h )*dΓ
u0(v) = ∫( v*u0f )*dΩ
massop = AffineFEOperator(m,u0,Ug,V0)
op = AffineFEOperator(a,b,Ug,V0)

Mass = get_matrix(massop)
y0 = get_vector(massop)
L = get_matrix(op)

function fx(t,ftrue,a,Ug,V0,h,dΩ,dΓ)
    f(x) = ftrue(x,t)
    b(v) = ∫( v*f )*dΩ + ∫( v*h )*dΓ
    op = AffineFEOperator(a,b,Ug,V0)
    return get_vector(op)
end

ff(t) = fx(t,ftrue,a,Ug,V0,h,dΩ,dΓ)

tspan = [0,1]
n = parse(Int64,ARGS[4])
println("n: ", n)

maxsylvit = 60
t,y,itersyl,ressyl = rk_lin_thr_solve_mex(method,s,Mass,L,ff,tspan,n,y0,maxsylvit)
@btime t,y,itersyl,ressyl = rk_lin_thr_solve_mex($method,$s,$Mass,$L,ff,$tspan,$n,$y0,$maxsylvit)

println("Sylvester convergence:")
for i in eachindex(itersyl)
    println("Iterations: ", itersyl[i], " Residual: ", ressyl[i])
end

uh = FEFunction(Ug,y[:,n])

# Create a string containing result_ followed by the method, s, n and the number of threads
filename = "result_" * method * "_" * string(s) * "_" * string(n) * "_" * string(nthreads()) 
writevtk(Ω,filename,cellfields=["uh"=>uh])
