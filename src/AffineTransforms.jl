__precompile__(true)

module AffineTransforms

import Base: *, \, eltype, ndims, show
if isdefined(Base, :scale)
    import Base.scale
else
    export scale
end

using Interpolations, Requires, Base.Cartesian

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
