using Base.Test

import AffineTransforms

# Pure translations
for nd = 1:4
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
end

# Pure rotations
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
# 2d
tfm1 = AffineTransforms.tformrotate(theta1)
tfm2 = AffineTransforms.tformrotate(theta1)
x = rand(2)
@test_approx_eq tfm1*(tfm2*x) -x
@test_approx_eq (tfm1*tfm2)*x -x
@test_approx_eq tfm1\(tfm2*x)  x
@test_approx_eq (tfm1\tfm2)*x  x
@test_approx_eq (tfm1*tfm1*tfm1*tfm1)*x  x

# 3d
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

# random transformations
for nd = 1:4
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
end

# Pure rotations in 2d and 3d
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
