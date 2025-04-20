"""
This function uses the IRK method to solve the ODE system
    M dy/dt = - L y + f(t)
    y(0) = y0
using shared memory parallelism 

Parameters
----------
method: String
    The name of the IRK method to use
s : Int
    The number of stages of the IRK method
Mass : Array{Float64, 2}
    The mass matrix of the ODE system
L : Array{Float64, 2}
    The stiffness matrix of the ODE system
f : Function
    The function f(t) that defines the ODE system
tspan : Tuple{Float64, Float64}
    The time interval [t0, tf] over which to solve the ODE system
n : Int
    The number of time steps to take
y0 : Array{Float64, 1}
    The initial condition of the ODE system
"""
function rk_lin_thr_solve_mex(method,s,Mass::SparseMatrixCSC,L::SparseMatrixCSC,f,tspan,n,y0,sylvester_maxit=50)

    # First we check that the number of spawned threads is less or equal than the number of stages of the method
    if nthreads() > s
        error("The number of threads must be less or equal to the number of stages of the method")
    end

    # First we generate the Butcher Tableau for the method
    # to employ the matrix equation solver inside the stage computation
    # we need the transpose of the matrix A for the Gauss scheme, and of the
    # the matrix M for the other schemes. We explicitly compute it here,
    # just once for all
    if method == "gauss"
        # For the Gauss method we use the fact that is is a symmetric scheme
        A, b, c, X, ev, D = generate_rkcoeff(s,method)
        At = convert(Matrix{Float64},copy(transpose(A)))
    else
        # In all the other cases we use the general form
        A, b, c, M, B, W, D = generate_rkcoeff(s,method)
        Mt = convert(Matrix{Float64},copy(transpose(M)))
    end 
    dss = D[s,s]
    
    # Build the time vector
    t0, tf = tspan
    h = (tf - t0)/n
    t = LinRange(t0,tf,n+1)
    # Allocate space for the solution matrix
    m = length(y0)
    y = zeros(m, n+1)
    y[:,1] = y0

    # Build a sparse identity matrix of size m
    es = ones(s,1)  

    # We will put here the evaluation of the function f at the current stage values in the Newton Iteration
    F = zeros(Float64,m,s);
    K = zeros(ComplexF64,m,s);
    E = zeros(Float64,m,s);         # Preallocated for speed
    Ksol = zeros(Float64,m,s);      # Preallocated for speed
    # Sylvester equation infos
    itersyl = zeros(Int64,n)
    ressyl = zeros(Float64,n)

    if method == "gauss"

        # Preallocate the matrices for the Sylvester equation
        VV = zeros(m,sylvester_maxit)
        HH = zeros(sylvester_maxit+1,sylvester_maxit)
        # Factorize matrices once for all
        Xt = copy(transpose(X))
        Xf = factorize(Xt)       # Factorize X
        Lfact = factorize(L)    # Factorize L 
        # Create a vector of length s with the factorization of Mass + h*ev[j]*L for j = 1:s 
        lhs = fetch.([Threads.@spawn factorize(Mass + h*ev[j]*L) for j in 1:s])
        vs = transpose(es)/Xf

        for i = 1:length(t)-1
            ti = t[i];
            # Assemble right-hand side
            Threads.@threads for j in eachindex(c)
                @inbounds F[:,j] = f(ti + c[j]*h);
            end
            
            r = F/Xf - (L*y[:,i])*vs;

            # Solve the block diagonal linear system
            Threads.@threads for j in eachindex(lhs)
                 @inbounds K[:,j] = lhs[j]\r[:,j]; 
            end
            
            # Compute the stages
            K *= Xt;
            
            # Convert K to Float64
            Ksol = real(K);

            # Solve the Sylvester equation
            E,itersyl[i],ressyl[i] = proj_sylvesterc_prealloc(h,Mass,Lfact,At,-Ksol*b/2.0,es,sylvester_maxit,1e-10,VV,HH)

            # Correct K 
            Ksol = Ksol + E;

            # Update the solution using an iterator to avoid race conditions
            temp = zeros(m)
            for j=1:s
                temp += h*b[j]*Ksol[:,j];
            end
            y[:,i+1] = y[:,i] + temp
        end

    elseif method == "radauIA" || method == "radauIIA"
        VV = zeros(m,sylvester_maxit*2)
        HH = zeros(sylvester_maxit*2+1,sylvester_maxit*2)
        # Here there is no correction on D

        C1 = zeros(s,2); C1[1,1] = -0.5; C1[s,2] = -1.0/(4.0*s - 2);
        C2 = zeros(s,2); C2[1,1] = 1.0; C2[s,2] = 1.0;

        Xhat = M + C1*transpose(C2)
        # Set the entries of Xhat that are small explicitly to zero
        Xhat = map(x -> abs(x) < 1e-13 ? 0 : x, Xhat) 
        ev, Q = eigen(Xhat)
        # Set the real part of ev to zeros
        ev = complex.(zeros(s), imag(ev))

        # Factorize matrices once for all
        # lhs = [factorize(Mass + h*ev[j]*L) for j=1:s]   # Factorize the blocks of the block diagonal matrix
        lhs = fetch.([Threads.@spawn factorize(Mass + h*ev[j]*L) for j in 1:s])
        Lfact = factorize(L)                            # Factorize L

        for i = 1:length(t)-1
            ti = t[i];
            # Assemble right-hand side
            Threads.@threads for j in eachindex(c)
                @inbounds F[:,j] = f(ti + c[j]*h);
            end

            r = (-L*y[:,i]*transpose(b)*W + F*B*W)*transpose(Q');

            # Solve the block diagonal linear system
            Threads.@threads for j in eachindex(lhs)
                 @inbounds K[:,j] = lhs[j]\r[:,j]; 
            end

            # Compute the stages
            K = K*transpose(Q);

            Ksol = real(K)

            # Solve the Sylvester equation
            E,itersyl[i],ressyl[i] = proj_sylvesterc_block_prealloc(h,Mass,Lfact,Mt,Ksol*C2,C1,m,1e-10,VV,HH)

            # Correct K 
            Ksol = Ksol + E;

            # Convert K to Float64 and go back to the original variables
            Ksol = real(Ksol)*transpose(W);    

            # Update the solution
            temp = zeros(m,1)
            for j=1:s
                temp = temp + h*b[j]*Ksol[:,j];
            end
            y[:,i+1] = y[:,i] + temp

        end    

    else
        # These methods have the diagonal D

        if method != "lobattoIIID"
            rkcor = 2
        else
            rkcor = 1
        end
        C1 = zeros(s,rkcor)
        C2 = zeros(s,rkcor)
        C1[1,1] = -0.5
        @inline xicoeff(k) = 1.0/(2.0*sqrt(4.0*k^2-1))
        if method == "lobattoIIIA"
            C1[s-1,2] = - ((2*s-1)/(s-1))*xicoeff(s-1)
            C2[1,1] = 1.0; C2[s,2] = 1.0;
        elseif method == "lobattoIIIB"
            C1[s,2] = ((2*s-1)/(s-1))*xicoeff(s-1)
            C2[1,1] = 1.0; C2[s-1,2] = 1.0;
        elseif method == "lobattoIIIC"
            C1[s,2] = - ((2.0*s-1.0)/((2.0*s-2)*(s-1)))
            C2[1,1] = 1.0; C2[s,2] = 1.0;
        elseif method == "lobattoIIID"
            C2[1,1] = 1.0
        end   
        
        # Preallocate the matrices for the Sylvester equation
        VV = zeros(m,sylvester_maxit*rkcor)
        HH = zeros(sylvester_maxit*rkcor+1,sylvester_maxit*rkcor)

        is = zeros(s); is[s] = 1.0;

        Xhat = M + C1*transpose(C2)
        # Set the entries of Xhat that are small explicitly to zero
        Xhat = map(x -> abs(x) < 1e-13 ? 0 : x, Xhat) 
        ev, Q = eigen(Xhat)
        # Set the real part of ev to zeros
        ev = complex.(zeros(s), imag(ev))

        # Factorize matrices once for all
        # lhs = [factorize(Mass + h*ev[j]*L) for j=1:s]   # Factorize the blocks of the block diagonal matrix
        lhs = fetch.([Threads.@spawn factorize(Mass + h*ev[j]*L) for j in 1:s])
        Lfact = factorize(L)                            # Factorize L

        for i = 1:length(t)-1
            ti = t[i];
            # Assemble right-hand side
            Threads.@threads for j in eachindex(c)
                @inbounds F[:,j] = f(ti + c[j]*h);
            end

            # r = reshape((-L*y[:,i]*transpose(b)*W + F*B*W)*transpose(Q'),m*s,1)
            r = (-L*y[:,i]*transpose(b)*W + F*B*W)*transpose(Q')

            # Solve the block diagonal linear system
            Threads.@threads for j in eachindex(lhs)
                 @inbounds K[:,j] = lhs[j]\r[:,j];
            end

            # Compute the stages
            K = K*transpose(Q);

            Ksol = real(K)

            # Build the right-hand side of the Sylvester equation
            cL = ((1-dss)/h)*(L\(Mass*(Ksol*is)))
            cR = D\is
            C2t = Ksol*C2
            C1t = D\C1

            C1hat,C2hat = build_rhs_qr(cL,cR,C2t,C1t)

            # Solve the Sylvester equation
            E,itersyl[i],ressyl[i] = proj_sylvesterc_block_prealloc(h,Mass,Lfact,Mt,C1hat,C2hat,m,1e-10,VV,HH)
            
            # Correct K 
            Ksol = Ksol + E;

            # Convert K to Float64 and go back to the original variables
            Ksol = real(Ksol)*transpose(W);    

            # Update the solution
            temp = zeros(m,1)
            for j=1:s
                temp = temp + h*b[j]*Ksol[:,j];
            end
            y[:,i+1] = y[:,i] + temp

        end    

    end

    return t,y,itersyl,ressyl

end

"""
This function uses the IRK method to solve the ODE system
    M dy/dt = G(t,y)
    y(0) = y0
using shared memory parallelism 

Parameters
----------
method: String
    The name of the IRK method to use
s : Int
    The number of stages of the IRK method
Mass : Array{Float64, 2}
    The mass matrix of the ODE system
Theta : Function
    The function Theta(t,y) that defines the ODE system
JTheta : Function
    The Jacobian of Theta(t,y)
tspan : Tuple{Float64, Float64}
    The time interval [t0, tf] over which to solve the ODE system
n : Int
    The number of time steps to take
y0 : Array{Float64, 1}
    The initial condition of the ODE system
"""
function rk_nlin_thr_solve_mex(method,s,Mass::SparseMatrixCSC,Theta,JTheta,tspan,n,y0,sylvester_maxit=50)

    # First we check that the number of spawned threads is less or equal than the number of stages of the method
    if nthreads() > s
        error("The number of threads must be less or equal to the number of stages of the method")
    end

    # First we generate the Butcher Tableau for the method
    # to employ the matrix equation solver inside the stage computation
    # we need the transpose of the matrix A for the Gauss scheme, and of the
    # the matrix M for the other schemes. We explicitly compute it here,
    # just once for all
    if method == "gauss"
        # For the Gauss method we use the fact that is is a symmetric scheme
        A, b, c, X, ev, D = generate_rkcoeff(s,method)
        At = convert(Matrix{Float64},copy(transpose(A)))
    else
        # In all the other cases we use the general form
        A, b, c, M, B, W, D = generate_rkcoeff(s,method)
        Mt = convert(Matrix{Float64},copy(transpose(M)))
    end 
    dss = D[s,s]
    
    # Build the time vector
    t0, tf = tspan
    h = (tf - t0)/n
    t = LinRange(t0,tf,n+1)
    # Allocate space for the solution matrix
    m = length(y0)
    y = zeros(m, n+1)
    y[:,1] = y0

    # Build a sparse identity matrix of size m
    es = ones(s,1)  

    # We will put here the evaluation of the function f at the current stage values in the Newton Iteration
    F = zeros(Float64,m,s);
    K = zeros(Float64,m,s);
    E = zeros(Float64,m,s);         # Preallocated for speed
    Dnewt = zeros(ComplexF64,m,s);  # Preallocated for speed
    Dnewt_sol = zeros(Float64,m,s); # Preallocated for speed
    iter = 0                        # Preallocated for speed
    res = 0.0                       # Preallocated for speed
    itersyl = zeros(Int64,n)
    ressyl = zeros(Float64,n)

    newton_maxit = 100
    tol_newton = 1e-10

    if method == "gauss"

        # Factorize matrices once for all
        Xt = copy(transpose(X))
        Xf = factorize(Xt)       # Factorize X
        Massfact = factorize(Mass) # Factorize Mass
        VV = zeros(m,sylvester_maxit*2)
        HH = zeros((sylvester_maxit+1)*2,sylvester_maxit*2)


        for i = 1:length(t)-1
            ti = t[i];

            # Use simplified Newton to solve for the stages
            L = JTheta(ti,y[:,i]) # Compute the Jacobian of Theta at the current time and state
            Lfact = factorize(L)  # Factorize L (I need it for the Sylvester equation)
            lhs = fetch.([Threads.@spawn factorize(Mass + h*ev[j]*L) for j in 1:s])
            K .= 0.0
            Am(x) = (Lfact\(Mass*x))
            As(x) = (Massfact\(L*x))
            println("Time step: ", i)
            for p = 1:newton_maxit
                # Assemble right-hand side
                for j in eachindex(c)
                    F[:,j] = Mass*(K[:,j] - y[:,i])  + h*sum(A[j,k]*Theta(ti + c[k]*h,K[:,k]) for k=1:s);
                end

                normF = norm(F)

                if  normF < tol_newton
                    println("Newton iteration: ", p, " Residual: ", normF)
                    itersyl[i] = round(itersyl[i]/p)
                    ressyl[i] = ressyl[i]/p
                    break
                end 

                r = F/Xf

                # Solve the Newton system
                Threads.@threads for j in eachindex(lhs)
                    @inbounds Dnewt[:,j] = lhs[j]\r[:,j]; 
                end

                # Compute the stages
                Dnewt = Dnewt*Xt;
                
                # Convert Dnewt_sol to Float64
                Dnewt_sol = real(Dnewt);

                # Solve the Sylvester equation
                E,iter,res = kpik_sylv_oneside_prealloc(Am, As, h*At, -h*Dnewt_sol*b/2.0,es, m, 1e-10,VV,HH)
                itersyl[i] = 2*itersyl[i] + iter
                ressyl[i] = ressyl[i] + res                 

                # Perform the Newton update on the K 
                K = K - (Dnewt_sol + E)

            end
            

            # Update the solution using an iterator to avoid race conditions
            temp = zeros(m)
            for j=1:s
                temp -= h*b[j]*Theta(ti + c[j]*h,K[:,j]);
            end
            y[:,i+1] = y[:,i] + temp
        end

    elseif method == "radauIA" || method == "radauIIA"
        # Here there is no correction on D

        C1 = zeros(s,2); C1[1,1] = -0.5; C1[s,2] = -1.0/(4.0*s - 2);
        C2 = zeros(s,2); C2[1,1] = 1.0; C2[s,2] = 1.0;
        VV = zeros(m,sylvester_maxit*4)
        HH = zeros((sylvester_maxit+1)*4,sylvester_maxit*4)

        Xhat = M + C1*transpose(C2)
        # Set the entries of Xhat that are small explicitly to zero
        Xhat = map(x -> abs(x) < 1e-13 ? 0 : x, Xhat) 
        ev, Q = eigen(Xhat)
        # Set the real part of ev to zeros
        ev = complex.(zeros(s), imag(ev))

        # Factorize matrices once for all
        lhs = [factorize(Mass + h*ev[j]*L) for j=1:s]   # Factorize the blocks of the block diagonal matrix
        Lfact = factorize(L)                            # Factorize L

        for i = 1:length(t)-1
            ti = t[i];
            # Assemble right-hand side
            for j=1:s
                F[:,j] = f(ti + c[j]*h);
            end

            r = reshape((-L*y[:,i]*transpose(b)*W + F*B*W)*transpose(Q'),m*s,1)

            # Solve the block diagonal linear system
            for j=1:s
                K[:,j] = lhs[j]\r[(j-1)*m+1:j*m]; 
            end

            # Compute the stages
            K = K*transpose(Q);

            Ksol = real(K)

            # Solve the Sylvester equation
            kpik_sylv_oneside_prealloc(Am, As, Mt, Ksol*C2, C1, m, 1e-10,VV,HH)

            # Correct K 
            Ksol = Ksol + E;

            # Convert K to Float64 and go back to the original variables
            Ksol = real(Ksol)*transpose(W);    

            # Update the solution
            temp = zeros(m,1)
            for j=1:s
                temp = temp + h*b[j]*Ksol[:,j];
            end
            y[:,i+1] = y[:,i] + temp

        end    

    else
        # These methods have the diagonal D

        if method != "lobattoIIID"
            rkcor = 2
        else
            rkcor = 1
        end
        C1 = zeros(s,rkcor)
        C2 = zeros(s,rkcor)
        C1[1,1] = -0.5
        @inline xicoeff(k) = 1.0/(2.0*sqrt(4.0*k^2-1))
        if method == "lobattoIIIA"
            C1[s-1,2] = - ((2*s-1)/(s-1))*xicoeff(s-1)
            C2[1,1] = 1.0; C2[s,2] = 1.0;
        elseif method == "lobattoIIIB"
            C1[s,2] = ((2*s-1)/(s-1))*xicoeff(s-1)
            C2[1,1] = 1.0; C2[s-1,2] = 1.0;
        elseif method == "lobattoIIIC"
            C1[s,2] = - ((2.0*s-1.0)/((2.0*s-2)*(s-1)))
            C2[1,1] = 1.0; C2[s,2] = 1.0;
        elseif method == "lobattoIIID"
            C2[1,1] = 1.0
        end   
        
        VV = zeros(m,sylvester_maxit*2*rkcor)
        HH = zeros((sylvester_maxit+1)*2*rkcor,sylvester_maxit*2*rkcor)

        is = zeros(s); is[s] = 1.0;

        Xhat = M + C1*transpose(C2)
        # Set the entries of Xhat that are small explicitly to zero
        Xhat = map(x -> abs(x) < 1e-13 ? 0 : x, Xhat) 
        ev, Q = eigen(Xhat)
        # Set the real part of ev to zeros
        ev = complex.(zeros(s), imag(ev))

        # Factorize matrices once for all
        lhs = [factorize(Mass + h*ev[j]*L) for j=1:s]   # Factorize the blocks of the block diagonal matrix
        Lfact = factorize(L)                            # Factorize L
        for i = 1:length(t)-1
            ti = t[i];
            # Assemble right-hand side
            for j=1:s
                F[:,j] = f(ti + c[j]*h);
            end

            r = reshape((-L*y[:,i]*transpose(b)*W + F*B*W)*transpose(Q'),m*s,1)

            # Solve the block diagonal linear system
            for j=1:s
                K[:,j] = lhs[j]\r[(j-1)*m+1:j*m]; 
            end

            # Compute the stages
            K = K*transpose(Q);

            Ksol = real(K)

            # Build the right-hand side of the Sylvester equation
            cL = ((1-dss)/h)*(L\(Mass*(Ksol*is)))
            cR = D\is
            C2t = Ksol*C2
            C1t = D\C1

            C1hat,C2hat = build_rhs_qr(cL,cR,C2t,C1t)

            # Solve the Sylvester equation
            E,iter,res = proj_sylvesterc_block(h,Mass,Lfact,Mt,C1hat,C2hat,m,1e-10)


            # Correct K 
            Ksol = Ksol + E;

            # Convert K to Float64 and go back to the original variables
            Ksol = real(Ksol)*transpose(W);    

            # Update the solution
            temp = zeros(m,1)
            for j=1:s
                temp = temp + h*b[j]*Ksol[:,j];
            end
            y[:,i+1] = y[:,i] + temp

        end    

    end

    return t,y,itersyl,ressyl

end