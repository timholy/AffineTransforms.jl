__precompile__(true)

module AffineTransforms

import Base: *, \, eltype, ndims, scale, show
using Interpolations, Requires, Base.Cartesian

if !isdefined(:AbstractVecOrMat)
    typealias AbstractVecOrMat Union(AbstractVector, AbstractMatrix)
end

export
    # types
    AffineTransform,
    TransformedArray,
    # functions
    rotation2,
    rotation3,
    rotationmatrix,
    rotationparameters,
    tformeye,
    tformfwd,
    tformfwd!,
    tforminv,
    tforminv!,
    tformfwd3rd,
    tformrigid,
    tformrotate,
    tformscale,
    tformtranslate,
    transform,
    transform!

include("transform.jl")
include("tformedarrays.jl")

end
