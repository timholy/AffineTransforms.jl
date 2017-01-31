using Base.Test

import AffineTransforms

# Pure translations
@testset "Pure $(nd)D translations" for nd = 1:4
    t1 = rand(nd)
    t2 = rand(nd)
    ident = eye(nd)

    tfm1 = AffineTransforms.AffineTransform(ident, t1)
    x = zeros(nd)
    @test AffineTransforms.tformfwd(tfm1, x) == t1
    @test_throws DimensionMismatch AffineTransforms.tformfwd(tfm1, zeros(nd+1))
    @test AffineTransforms.tforminv(tfm1, x) == -t1
    x = rand(nd)
    @test AffineTransforms.tformfwd(tfm1, x) == x+t1
    @test AffineTransforms.tforminv(tfm1, x) == x-t1
    @test tfm1*x == x+t1
    @test tfm1\x == x-t1
    X = reshape(x, nd, 1)
    @test all(AffineTransforms.tformfwd(tfm1, X) .== x+t1)
    @test all(AffineTransforms.tforminv(tfm1, X) .== x-t1)

    tfm2 = AffineTransforms.tformtranslate(t2)
    tfmP1 = tfm1*tfm2
    @test tfmP1.offset == t1+t2
    @test_approx_eq AffineTransforms.tformfwd(tfmP1, x) x+t2+t1

    if 2 <= nd <= 3
        # Tuple interface
        r = AffineTransforms.tformfwd(tfm1, x...)
        @test [r...] == tfm1*x
        r = AffineTransforms.tforminv(tfm1, x...)
        @test [r...] == tfm1\x
    end
end

@testset "Pure rotations" begin
    @test_approx_eq AffineTransforms.tformrotate(0).scalefwd eye(2)
    @test_approx_eq AffineTransforms.tformrotate(pi/2)*[1,0] [0,1]
    @test_approx_eq AffineTransforms.tformrotate(pi/2).scalefwd [0 -1; 1 0]
    @test_approx_eq AffineTransforms.tformrotate(pi)*[1,0] [-1,0]
    @test_approx_eq AffineTransforms.tformrotate(3pi/2)*[1,0] [0,-1]
    @test_approx_eq AffineTransforms.tformrotate([1,0,0],pi/2).scalefwd [1 0 0; 0 0 -1; 0 1 0]
    @test_approx_eq AffineTransforms.tformrotate([1,0,0],pi/2)*[1,0,0] [1,0,0]
    @test_approx_eq AffineTransforms.tformrotate([0,1,0],pi/2)*[1,0,0] [0,0,-1]
    @test_approx_eq AffineTransforms.tformrotate([0,0,1],pi/2)*[1,0,0] [0,1,0]
    s2 = 1/sqrt(2)
    @test_approx_eq AffineTransforms.tformrotate([0,1,0],pi/4).scalefwd [s2 0 s2; 0 1 0; -s2 0 s2]
    @test_approx_eq AffineTransforms.tformrotate([0,0,1],pi/4).scalefwd [s2 -s2 0; s2 s2 0; 0 0 1]

    theta1 = pi/2
    @testset "> 2d" begin
        tfm1 = AffineTransforms.tformrotate(theta1)
        tfm2 = AffineTransforms.tformrotate(theta1)
        x = rand(2)
        @test_approx_eq tfm1*(tfm2*x) -x
        @test_approx_eq (tfm1*tfm2)*x -x
        @test_approx_eq tfm1\(tfm2*x)  x
        @test_approx_eq (tfm1\tfm2)*x  x
        @test_approx_eq (tfm1*tfm1*tfm1*tfm1)*x  x
        r = AffineTransforms.tformfwd(tfm1, x...)
        @test [r...] == tfm1*x
        r = AffineTransforms.tforminv(tfm1, x...)
        @test [r...] == tfm1\x
    end

    @testset "> 3d" begin
        uaxis = rand(3); uaxis = uaxis/norm(uaxis)
        tfm1 = AffineTransforms.tformrotate(uaxis*theta1)
        tfm2 = AffineTransforms.tformrotate(uaxis*theta1)
        x = rand(3)
        xparallel = dot(x,uaxis)*uaxis
        xperp = x - xparallel
        @test_approx_eq tfm1*(tfm2*x) xparallel-xperp
        @test_approx_eq (tfm1*tfm2)*x xparallel-xperp
        @test_approx_eq tfm1\(tfm2*x)  x
        @test_approx_eq (tfm1\tfm2)*x  x
        @test_approx_eq (tfm1*tfm1*tfm1*tfm1)*x  x
        r = AffineTransforms.tformfwd(tfm1, x...)
        @test [r...] == tfm1*x
        r = AffineTransforms.tforminv(tfm1, x...)
        @test [r...] == tfm1\x
    end
end

@testset "Random $(nd)D transformations" for nd = 1:4
    A = eye(nd) + 0.1*randn(nd,nd)
    b = randn(nd)
    tfm = AffineTransforms.AffineTransform(A, b)
    x = rand(nd)
    @test_approx_eq tfm*x A*x+b
    @test_approx_eq tfm\x A\(x-b)
    @test_approx_eq (tfm\tfm)*x x
    # Scaling transformations
    s = rand(nd)
    y = tfm*x
    @test_approx_eq scale(s, tfm)*x  s.*y
    @test_approx_eq scale(tfm, s)*x  tfm*(s.*x)
    if 2 <= nd <= 3
        r = AffineTransforms.tformfwd(tfm, x...)
        @test [r...] == tfm*x
        r = AffineTransforms.tforminv(tfm, x...)
        @test [r...] == tfm\x
    end
end

@testset "Pure rotations in 2d and 3d" begin
    for i = 1:20
        th = pi*(rand()-0.5)
        R = AffineTransforms.rotation2(th)
        @test_approx_eq th AffineTransforms.rotationparameters(R)
    end
    for i = 1:100
        ax = randn(3)
        n = norm(ax)
        if n == 0
            continue
        end
        ax = ax / n
        th = pi*rand()
        rv = th*ax
        R = AffineTransforms.rotation3(rv)
        @test_approx_eq rv AffineTransforms.rotationparameters(R)
    end
end

@testset "Transformation of arrays" begin
    using Interpolations, Images # Images for padarray
    import AffineTransforms: center

    padwith0(A) = parent(padarray(A, Fill(0, (2,2))))

    for (IT,GT) in ((BSpline(Constant()), OnCell()),
                    (BSpline(Linear()), OnGrid()),
                    (BSpline(Quadratic(Flat())), OnCell()))
        A = padwith0([1 1; 1 2])
        itp = interpolate(A, IT, GT)
        tfm = AffineTransforms.tformtranslate([1,0])
        tA = AffineTransforms.TransformedArray(extrapolate(itp, NaN), tfm)
        @test_approx_eq AffineTransforms.transform(tA)[3:4,3:4] [1 2; 0 0]
        @test_approx_eq tA[3,3] 1
        @test_approx_eq tA[3,4] 2
        tfm = AffineTransforms.tformtranslate([0,1])
        tA = AffineTransforms.TransformedArray(extrapolate(itp, NaN), tfm)
        @test_approx_eq AffineTransforms.transform(tA)[3:4,3:4] [1 0; 2 0]
        tfm = AffineTransforms.tformtranslate([-1,0])
        tA = AffineTransforms.TransformedArray(extrapolate(itp, NaN), tfm)
        @test_approx_eq AffineTransforms.transform(tA)[3:4,3:4] [0 0; 1 1]
        tfm = AffineTransforms.tformtranslate([0,-1])
        tA = AffineTransforms.TransformedArray(extrapolate(itp, NaN), tfm)
        @test_approx_eq AffineTransforms.transform(tA)[3:4,3:4] [0 1; 0 1]

        A = padwith0([2 1; 1 1])
        itp = interpolate(A, IT, GT)
        tfm = AffineTransforms.tformrotate(pi/2)
        tA = AffineTransforms.TransformedArray(extrapolate(itp, NaN), tfm)
        @test_approx_eq AffineTransforms.transform(tA)[3:4,3:4] [1 2; 1 1]
        tfm = AffineTransforms.tformrotate(-pi/2)
        tA = AffineTransforms.TransformedArray(extrapolate(itp, NaN), tfm)
        @test_approx_eq AffineTransforms.transform(tA)[3:4,3:4] [1 1; 2 1]

        # Check that getindex and transform yield the same results
        A = rand(6,6)
        itp = interpolate(A, IT, GT)
        tfm = AffineTransforms.tformrotate(pi/7)*AffineTransforms.tformtranslate(rand(2))
        tA = AffineTransforms.TransformedArray(extrapolate(itp, NaN), tfm)
        dest = AffineTransforms.transform(tA)
        tfm_recentered = AffineTransforms.AffineTransform(tfm.scalefwd, tfm.offset + center(A) - tfm.scalefwd*center(dest))
        tA_recentered = AffineTransforms.TransformedArray(extrapolate(itp, NaN), tfm_recentered)
        nbad = 0
        for j = 1:size(dest,2), i = 1:size(dest,1)
            val = tA_recentered[i,j]
            # try
            #     @test isfinite(dest[i,j]) == isfinite(val)
            # catch err
            #     @show dest tA_recentered
            #     rethrow(err)
            # end
            nbad += isfinite(dest[i,j]) != isfinite(val)
            if isfinite(val) && isfinite(dest[i,j])
                @test abs(val-dest[i,j]) < 1e-12
            end
        end
        @test nbad < 3
        dest = AffineTransforms.transform(tA, origin_src=zeros(2), origin_dest=zeros(2))
        nbad = 0
        for j = 1:size(dest,2), i = 1:size(dest,1)
            val = tA[i,j]
            # try
            #     @test isfinite(dest[i,j]) == isfinite(val)
            # catch err
            #     @show dest tA
            #     rethrow(err)
            # end
            nbad += isfinite(dest[i,j]) != isfinite(val)
            if isfinite(val) && isfinite(dest[i,j])
                @test abs(val-dest[i,j]) < 1e-12
            end
        end
        @test nbad < 3
    end

    A = Float64[1 2 3 4; 5 6 7 8]
    tfm = AffineTransforms.tformeye(2)
    dest = zeros(2,2)

    itp = interpolate(A, BSpline(Linear()), OnGrid())
    tA = AffineTransforms.TransformedArray(extrapolate(itp, NaN), tfm)
    @test_approx_eq AffineTransforms.transform!(dest, tA) [2 3; 6 7]
    dest3 = zeros(3,3)
    @test_approx_eq AffineTransforms.transform!(dest3, tA) [NaN NaN NaN; 3.5 4.5 5.5; NaN NaN NaN]

    itp = interpolate(A, BSpline(Constant()), OnCell())
    tA = AffineTransforms.TransformedArray(extrapolate(itp, NaN), tfm)
    @test AffineTransforms.transform!(dest, tA) == [2 3; 6 7]

    # Transforms with non-real numbers
    using DualNumbers
    a = AffineTransforms.tformtranslate([0,dual(1,0)])
    A = reshape(1:9, 3, 3)
    At = AffineTransforms.transform(A, a)
    @test At[:,1:2] == A[:,2:3]
end
