VERSION >= v"0.4.0-dev+6521" && __precompile__(true)

module AffineTransforms

import Base: *, \, eltype, ndims, scale, show

if !isdefined(:AbstractVecOrMat)
    typealias AbstractVecOrMat Union(AbstractVector, AbstractMatrix)
end

export
    # types
    AffineTransform,
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
    tformtranslate

include("transform.jl")

end
