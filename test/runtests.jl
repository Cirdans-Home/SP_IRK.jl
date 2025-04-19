using SP_IRK
using Test

@testset "SP_IRK.jl" begin
    
    # Test the generation of RK coefficients via the generate_rkcoeff function
    @testset "generate_rkcoeff" begin
        # Test for Gauss-Legendre
        A, b, c, X, ev, D = generate_rkcoeff(3, "gauss")
        @test size(A) == (3, 3)
        @test size(b) == (3,)
        @test size(c) == (3,)
        @test size(X) == (3, 3)
        @test size(ev) == (3,)
        @test size(D) == (3, 3)

        # Test for RadauIA
        A, b, c, M, B, W, D = generate_rkcoeff(3, "radauIA")
        @test size(A) == (3, 3)
        @test size(b) == (3,)
        @test size(c) == (3,)
        @test size(M) == (3, 3)
        @test size(B) == (3, 3)
        @test size(W) == (3, 3)
        @test size(D) == (3, 3)

        # Test for RadauIIA
        A, b, c, M, B, W, D = generate_rkcoeff(4, "radauIIA")
        @test size(A) == (4, 4)
        @test size(b) == (4,)
        @test size(c) == (4,)
        @test size(M) == (4, 4)
        @test size(B) == (4, 4)
        @test size(W) == (4, 4)
        @test size(D) == (4, 4)

    end

end
