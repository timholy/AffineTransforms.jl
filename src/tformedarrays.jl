import Interpolations: AbstractExtrapolation

"""
A `TransformedArray` is an `AbstractArray` type with an affine
coordinate shift: `A[i,j]` evaluates the "parent" array `P` used to
construct `A` at a location `x,y` given by `[x, y] = tfm*[i,j]`.  It
therefore allows you to use lazy-evaluation to perform affine
coordinate transformations.

```
A = TransformedArray(etp, tfm)
```
where `etp` is an extrapolation object as defined by the
Interpolations package, and `tfm` is an `AffineTransform`.
"""
immutable TransformedArray{T,N,E<:AbstractExtrapolation,Tf<:AffineTransform} <: AbstractArray{T,N}
    data::E
    tform::Tf
end
TransformedArray{T,N}(etp::AbstractExtrapolation{T,N}, a::AffineTransform) =
    TransformedArray{T,N,typeof(etp),typeof(a)}(etp, a)

function TransformedArray(A::AbstractInterpolation, a::AffineTransform)
    etp = extrapolate(A, NaN)
    TransformedArray(etp, a)
end

function TransformedArray(A::AbstractArray, a::AffineTransform)
    itp = interpolate(A, BSpline(Linear()), OnGrid())
    TransformedArray(itp, a)
end


Base.size(A::TransformedArray) = size(A.data)

function Base.getindex{T}(A::TransformedArray{T,2}, i::Number, j::Number)
    x, y = tformfwd(A.tform, i, j)
    A.data[x, y]
end

function Base.getindex{T}(A::TransformedArray{T,3}, i::Number, j::Number, k::Number)
    x, y, z = tformfwd(A.tform, i, j, k)
    A.data[x, y, z]
end

Base.similar{T}(A::TransformedArray, ::Type{T}, dims::Dims) = Array(T, dims)

"""
`transform(A, tfm; origin_dest=center(A), origin_src=center(A)`
computes the transformed `A` over its entire domain.  By default the
transformation is assumed to operate around the center of the input
array, and output coordinates are referenced relative to the center of
the output.

If `A` is a TransformedArray, then the syntax is just `transform(A;
origin_dest=center(A), origin_src=center(A))`. This is different from
the behavior of `A[:,:]`, which assumes the origin of coordinates to
be all-zeros.  To obtain behavior equivalent to `getindex`, supply
zero-vectors for both of them. Alternatively to make `getindex` behave
as `transform`, offset the origin of the transform used to construct
`A` by
```
origin_src - tform.scalefwd*origin_dest
```
"""
function transform(A::TransformedArray; kwargs...)
    dest = Array(promote_type(eltype(A.data), eltype(A.tform)), size(A))
    transform!(dest, A; kwargs...)
    dest
end

transform(A, a::AffineTransform; kwargs...) = transform(TransformedArray(A, a); kwargs...)

"""
`transform!(dest, src, tfm; origin_dest=center(dest),
origin_src=center(src))` is like `transform`, but using a
pre-allocated output array `dest`.

If `src` is already a TransformedArray, use `transform!(dest, src;
kwargs...)`.
"""
function transform!{S,T,N}(dest::AbstractArray{S,N},
                           src::TransformedArray{T,N};
                           origin_dest = center(dest),
                           origin_src = center(src))
    tform = src.tform
    if tform.scalefwd == eye(T, N) && tform.offset == zeros(N) && size(dest) == size(src) && origin_dest == origin_src
        copy!(dest, src)
        return dest
    end
    offset = tform.offset - tform.scalefwd*origin_dest + origin_src
    _transform!(dest, src, offset)
end

transform!(dest, src, a::AffineTransform; kwargs...) = transform!(dest, TransformedArray(src, a); kwargs...)

@require Images begin
    transform(A::Images.AbstractImage, a::AffineTransform; kwargs...) = Images.copyproperties(A, transform(Images.data(A), a; kwargs...))
    transform!(dest, A::Images.AbstractImage, a::AffineTransform; kwargs...) = Images.copyproperties(A, transform!(dest, Images.data(A), a; kwargs...))
end

# For a FilledExtrapolation, this is designed to (usually) avoid evaluating
# the interpolation unless it is in-bounds.  This often improves performance.
@generated function _transform!{S,T,N,E<:Interpolations.FilledExtrapolation}(
                                dest::AbstractArray{S,N},
                                src::TransformedArray{T,N,E},
                                offset)
    # Initialize the final column of s "matrix," e.g., s_1_3 = A_1_3 + o_3
    # s stands for source-coordinates. The first "column," s_d_1, corresponds
    # to the actual interpolation position. The later columns simply cache
    # previous computations.
    sN = ntuple(i->Expr(:(=), symbol(string("s_",i,"_$N")), Expr(:call, :+, symbol(string("A_",i,"_$N")), symbol(string("o_",i)))), N)
    quote
        tform = src.tform
        data = src.data
        @nexprs $N d->(o_d = offset[d])
        @nexprs $N j->(@nexprs $N i->(A_i_j = tform.scalefwd[i,j]))
        fill!(dest, data.fillvalue)
        $(sN...)
        @nloops($N,i,d->(d>1 ? (1:size(dest,d)) : (imin:imax)),
                # The pre-expression chooses the range within each column that will be in-bounds
                d->(d >  2 ? (@nexprs $N e->s_e_{d-1}=A_e_{d-1}+s_e_d) :
                    d == 2 ? begin
                        imin = 1
                        imax = size(dest,1)
                        @nexprs $N e->((imin, imax) = irange(imin, imax, A_e_1, s_e_2, size(data, e)))
                        @nexprs $N e->(s_e_1=A_e_1*imin+s_e_d)
                    end :
                    nothing), # pre
            d->(@nexprs $N e->(s_e_d += A_e_d)), # post
            # Perform the interpolation of the source data
            @inbounds (@nref $N dest i) = (@nref $N data d->s_d_1)
        )
        dest
    end
end

center(A::AbstractArray) = [(size(A,d)+1)/2 for d = 1:ndims(A)]

# Find i such that
#      imin <= i <= imax
#      1 <= coef*i+offset <= upper
function irange(imin::Int, imax::Int, coef, offset, upper)
    const thresh = 10^4/typemax(Int)  # needed to avoid InexactError for results with abs() bigger than typemax
    if real(coef) > thresh
        return max(imin, floor(Int, real((1-offset)/coef))), min(imax, ceil(Int, real((upper-offset)/coef)))
    elseif real(coef) < -thresh
        return max(imin, floor(Int, real((upper-offset)/coef))), min(imax, ceil(Int, real((1-offset)/coef)))
    else
        if 1 <= real(offset) <= upper
            return imin, imax
        else
            return 1, 0   # empty range
        end
    end
end
