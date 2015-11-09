immutable AffineTransform{T,N}
    scalefwd::Matrix{T}
    scaleinv::Matrix{T}
    offset::Vector{T}
    temp::Vector{T}

    function AffineTransform(S::AbstractMatrix, t::AbstractVector)
        if size(S) != (N,N) || size(t) != (N,)
            throw(DimensionMismatch("Inputs must be $N-dimensional"))
        end
        temp = Array(eltype(t), length(t))
        new(S, inv(S), t, temp)
    end
end
function AffineTransform(S, offset)
    T = promote_type(eltype(S),eltype(offset))
    AffineTransform{T,length(offset)}(convert(Matrix{T},S), convert(Vector{T},offset))
end

ndims{T,N}(tform::AffineTransform{T,N}) = N
eltype{T}(tf::AffineTransform{T}) = T

function show(io::IO, tf::AffineTransform)
    println(io, typeof(tf), ":")
    println(io, "matrix: ", tf.scalefwd)
    println(io, "translation: ", tf.offset)
end

scale(tf::AffineTransform, s::Vector) = AffineTransform(scale(tf.scalefwd, s), tf.offset)
scale(s::Vector, tf::AffineTransform) = AffineTransform(scale(s, tf.scalefwd), tf.offset.*s)

*(a::AffineTransform, v::AbstractVector) = tformfwd(a, v)
\(a::AffineTransform, x::AbstractVector) = tforminv(a, x)

tformfwd{T}(a::AffineTransform{T}, x::AbstractVector) = tformfwd!(Array(typeof(one(T)*one(eltype(x))), length(x)), a, x)
tforminv{T}(a::AffineTransform{T}, x::AbstractVector) = tforminv!(Array(typeof(one(T)*one(eltype(x))), length(x)), a, x)
tformfwd{T}(a::AffineTransform{T}, X::AbstractMatrix) = tformfwd!(Array(typeof(one(T)*one(eltype(X))), size(X,1), size(X,2)), a, X)
tforminv{T}(a::AffineTransform{T}, X::AbstractMatrix) = tforminv!(Array(typeof(one(T)*one(eltype(X))), size(X,1), size(X,2)), a, X)

tformfwd(a::AffineTransform, v::AbstractMatrix) = a.scalefwd*v .+ a.offset
tforminv(a::AffineTransform, x::AbstractMatrix) = a.scaleinv*(x .- a.offset)

# Implementation of a.scalefwd*X .+ a.offset
function tformfwd!(dest, a::AffineTransform, X::AbstractVecOrMat)
    A_mul_B!(dest, a.scalefwd, X)
    @inbounds for j = 1:size(X,2), i = 1:size(X,1)
        dest[i,j] += a.offset[i]
    end
    dest
end

# Implementation of a.scaleinv*(X .- a.offset)
function tforminv!(dest, a::AffineTransform, X::AbstractVecOrMat)
    size(X,1) == ndims(a) || throw(DimensionMismatch("Dimensionality $(size(X,1)) of inputs does not match $(ndims(a)) for the transformation"))
    A_mul_B!(dest, a.scaleinv, X)
    A_mul_B!(a.temp, a.scaleinv, a.offset)
    @inbounds for j = 1:size(X,2), i = 1:size(X,1)
        dest[i,j] -= a.temp[i]
    end
    dest
end

# Faster versions for particular dimensionalities (with vectors)
for VecType in (AbstractVector, CartesianIndex)
    @eval begin
        function tforminv!{T}(v, a::AffineTransform{T,3}, x::$VecType)
            length(x) == 3 || throw(DimensionMismatch("vector must be 3-dimensional"))
            t = a.temp
            o = a.offset
            @inbounds t[1] = x[1]-o[1]
            @inbounds t[2] = x[2]-o[2]
            @inbounds t[3] = x[3]-o[3]
            gemv3!(v, a.scaleinv, t)
        end
        function tformfwd!{T}(x, a::AffineTransform{T,3}, v::$VecType)
            length(x) == 3 || throw(DimensionMismatch("vector must be 3-dimensional"))
            gemv3!(x, a.scalefwd, v)
            o = a.offset
            @inbounds x[1] += o[1]
            @inbounds x[2] += o[2]
            @inbounds x[3] += o[3]
            x
        end
    end
end

function tforminv{T}(a::AffineTransform{T,3}, x1::Number, x2::Number, x3::Number)
    o = a.offset
    @inbounds t1 = x1 - o[1]
    @inbounds t2 = x2 - o[2]
    @inbounds t3 = x3 - o[3]
    gemv3(a.scaleinv, t1, t2, t3)
end

function tformfwd{T}(a::AffineTransform{T,3}, x1::Number, x2::Number, x3::Number)
    t1, t2, t3 = gemv3(a.scalefwd, x1, x2, x3)
    o = a.offset
    @inbounds r1 = t1+o[1]
    @inbounds r2 = t2+o[2]
    @inbounds r3 = t3+o[3]
    r1, r2, r3
end

for VecType in (AbstractVector, CartesianIndex)
    @eval begin
        function tforminv!{T}(v, a::AffineTransform{T,2}, x::$VecType)
            length(x) == 2 || throw(DimensionMismatch("vector must be 2-dimensional"))
            t = a.temp
            o = a.offset
            @inbounds t[1] = x[1]-o[1]
            @inbounds t[2] = x[2]-o[2]
            gemv2!(v, a.scaleinv, t)
        end
        function tformfwd!{T}(x, a::AffineTransform{T,2}, v::$VecType)
            length(x) == 2 || throw(DimensionMismatch("vector must be 2-dimensional"))
            gemv2!(x, a.scalefwd, v)
            t = a.offset
            @inbounds x[1] += t[1]
            @inbounds x[2] += t[2]
            x
        end
    end
end

function tforminv{T}(a::AffineTransform{T,2}, x1::Number, x2::Number)
    o = a.offset
    @inbounds t1 = x1 - o[1]
    @inbounds t2 = x2 - o[2]
    gemv2(a.scaleinv, t1, t2)
end

function tformfwd{T}(a::AffineTransform{T,2}, x1::Number, x2::Number)
    t1, t2 = gemv2(a.scalefwd, x1, x2)
    o = a.offset
    @inbounds r1 = t1+o[1]
    @inbounds r2 = t2+o[2]
    r1, r2
end

# For just the 3rd coordinate (we sometimes only need that one)
function tforminv3rd{T}(a::AffineTransform{T,3}, x::AbstractVector)
    o = a.offset
    A = a.scaleinv
    A[3]*(x[1]-o[1]) + A[6]*(x[2]-o[2]) + A[9]*(x[3]-o[3])
end

function scaled_add3!(x, s::Number, r)
    @inbounds x[1]+=s*r[1]
    @inbounds x[2]+=s*r[2]
    @inbounds x[3]+=s*r[3]
end

function gemv3!(res, A, x)
    @inbounds res[1] = A[1]*x[1] + A[4]*x[2] + A[7]*x[3]
    @inbounds res[2] = A[2]*x[1] + A[5]*x[2] + A[8]*x[3]
    @inbounds res[3] = A[3]*x[1] + A[6]*x[2] + A[9]*x[3]
    res
end

function gemv3(A, x1::Number, x2::Number, x3::Number)
    @inbounds r1 = A[1]*x1 + A[4]*x2 + A[7]*x3
    @inbounds r2 = A[2]*x1 + A[5]*x2 + A[8]*x3
    @inbounds r3 = A[3]*x1 + A[6]*x2 + A[9]*x3
    r1, r2, r3
end

function scaled_add2!(x, s::Number, r)
    @inbounds x[1]+=s*r[1]
    @inbounds x[2]+=s*r[2]
end

function gemv2!(res, A, x)
    @inbounds res[1] = A[1]*x[1] + A[3]*x[2]
    @inbounds res[2] = A[2]*x[1] + A[4]*x[2]
    res
end

function gemv2(A, x1::Number, x2::Number)
    @inbounds r1 = A[1]*x1 + A[3]*x2
    @inbounds r2 = A[2]*x1 + A[4]*x2
    r1, r2
end

### Identity transformation
tformeye(T, nd) = AffineTransform(eye(T, nd), zeros(T, nd))
tformeye(nd) = AffineTransform(eye(nd), zeros(nd))

### Rotations
# Right-hand rotation around an axis (axis points towards viewer)

rotation2(angle) = [cos(angle) -sin(angle); sin(angle) cos(angle)]

function tformrotate(angle)
    A = rotation2(angle)
    AffineTransform(A, zeros(eltype(A),2))
end

# The following assumes uaxis is normalized
function _rotation3(uaxis::Vector, angle)
    if length(uaxis) != 3
        error("3d rotations only")
    end
    ux, uy, uz = uaxis[1], uaxis[2], uaxis[3]
    c = cos(angle)
    s = sin(angle)
    cm = one(typeof(c)) - c
    R = [c+ux*ux*cm       ux*uy*cm-uz*s      ux*uz*cm+uy*s;
         uy*ux*cm+uz*s      c+uy*uy*cm       uy*uz*cm-ux*s;
         uz*ux*cm-uy*s    uz*uy*cm+ux*s      c+uz*uz*cm]
end

function rotation3{T}(axis::Vector{T}, angle)
    n = norm(axis)
    axisn = n>0 ? axis/n : (tmp = zeros(T,length(axis)); tmp[1] = 1)
    _rotation3(axisn, angle)
end

# Angle/axis representation where the angle is the norm of the vector (so axis is not normalized)
function rotation3{T}(axis::Vector{T})
    n = norm(axis)
    axisn = n>0 ? axis/n : (tmp = zeros(typeof(one(T)/1),length(axis)); tmp[1] = 1; tmp)
    _rotation3(axisn, n)
end

# Convert rotation matrix to axis representation (angle/axis with the norm of the vector coding the angle)
function rotationparameters{T}(R::Matrix{T})
    size(R, 1) == size(R, 2) || error("Matrix must be square")
    if size(R, 1) == 2
        return [atan2(-R[1,2],R[1,1])]
    elseif size(R, 1) == 3
        # See Itzhack Y. Bar-Itzhack.  "New Method for Extracting the Quaternion from a Rotation Matrix",
        # Journal of Guidance, Control, and Dynamics, Vol. 23, No. 6 (2000), pp. 1085-1087.
        # On https://en.wikipedia.org/wiki/Rotation_matrix#Conversions
        # This is said to be robust even if R is not orthogonal
        Q12p, Q13p, Q23p = R[1,2]+R[2,1], R[1,3]+R[3,1], R[2,3]+R[3,2]
        Q12m, Q13m, Q23m = R[1,2]-R[2,1], R[1,3]-R[3,1], R[2,3]-R[3,2]
        K = [R[1,1]-R[2,2]-R[3,3]  Q12p  Q13p  Q23m;
             Q12p  R[2,2]-R[1,1]-R[3,3]  Q23p -Q13m;
             Q13p  Q23p  R[3,3]-R[1,1]-R[2,2]  Q12m;
             Q23m -Q13m  Q12m  R[1,1]+R[2,2]+R[3,3]]
        D, V = eig(K)
        imax = indmax(D)
        q = V[:,imax]
        ax = q[1:3]
        n = norm(ax)
        if n > 0
            ax = ax/n
        end
        th = -2*atan(n/q[4])
        return th*ax
    else
        error("Rotations in $(size(R, 1)) dimensions not supported")
    end
end


tformtranslate(trans::Vector) = AffineTransform(eye(length(trans)), trans)

tformscale(scale::Real, ndims::Int) = AffineTransform(scale*eye(typeof(scale),ndims), zeros(typeof(scale),ndims))

function tformrotate(axis::Vector, angle)
    if length(axis) == 3
        return AffineTransform(rotation3(axis, angle), zeros(eltype(axis),3))
    else
        error("Dimensionality ", length(axis), " not supported")
    end
end

function tformrotate(x::Vector)
    if length(x) == 3
        return AffineTransform(rotation3(x), zeros(eltype(x),3))
    else
        error("Dimensionality ", length(x), " not supported")
    end
end

# The following is useful in optimization routines, which tend to work on parameter vectors
function tformrigid(p::Vector)
    if length(p) == 3
        # 2d: [angle, shiftx, shifty]
        return AffineTransform(rotation2(p[1]), p[2:3])
    elseif length(p) == 6
        # 3d: [axis, shift]
        return AffineTransform(rotation3(p[1:3]), p[4:6])
    else
        error("Rigid transformations with $(length(p)) parameters not supported")
    end
end

function(*)(T2::AffineTransform, T1::AffineTransform)
    A1 = T1.scalefwd
    A2 = T2.scalefwd
    o = T2.offset + A2*T1.offset
    AffineTransform(A2*A1, o)
end

function(\)(T2::AffineTransform, T1::AffineTransform)
    A1 = T1.scalefwd
    A2 = T2.scaleinv
    o = A2*(T1.offset-T2.offset)
    AffineTransform(A2*A1, o)
end
