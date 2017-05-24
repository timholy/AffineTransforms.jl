# AffineTransforms

[![Build Status](https://travis-ci.org/timholy/AffineTransforms.jl.svg?branch=master)](https://travis-ci.org/timholy/AffineTransforms.jl)

A package for working with affine transformations. For new projects, I recommend [CoordinateTransformations](https://github.com/FugroRoames/CoordinateTransformations.jl) instead.

## Installation

In julia, type
```julia
Pkg.add("AffineTransforms")
```

## Theory

An affine transformation is of the form
```julia
y = A*x + b
```
This is the "forward" transformation. The "inverse" transformation is therefore
```julia
x = A\(y-b)
```

## Usage

Create an affine transformation with
```julia
tfm = AffineTransform(A, b)
```
The following are all different ways of computing the forward transform:
```julia
y = tfm * x
y = tformfwd(tfm, x)
y = similar(x); tformfwd!(y, tfm, x)
```
Similarly, the following are all different ways of computing the inverse transform:
```julia
x = tfm\y
x = tforminv(tfm, y)
x = similar(y); tforminv!(x, tfm, y)
```

### Convenience constructors

```julia
tformeye(T, nd)
tformeye(nd)
```
Creates the identity transformation in `nd` dimensions.

```julia
tformtranslate(v)
```
Creates a shift (translation) transformation

```julia
tformrotate(angle)   # creates a 2d rotation
tformrotate(axis, angle)   # creates a 3d rotation
tformrotate(axis)          # creates a 3d rotation
```
In 3d, these constructors work with angle-axis representation, where `axis` is a 3-vector.
When `angle` is provided, `axis` is used as if it were normalized to have unit length.
If you just specify `axis`, then `norm(axis)` is used for the `angle`.

```julia
tformscale(scale::Real, nd)
```
Creates a scaling transformation, where `A` will have `scale` along the diagonal.

```julia
tformrigid(p)
```
Particularly useful for optimization of rigid transformations.
If `length(p) == 3`, this creates a 2d transform, where `p[1]` is the rotation angle, `p[2:3]` are
the two components of translation.
If `length(p) == 6`, this creates a 3d transform, where `p[1:3]` is `axis` for `tformrotate`,
and `p[4:6]` are the three components of translation.


### Representation conversions

```julia
rotationparameters(R)
```
Converts a 2d or 3d rotation matrix `R` into an `angle` (in 2d) or the `axis` representation (in 3d).
