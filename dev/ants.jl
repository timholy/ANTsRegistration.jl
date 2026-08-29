using Pkg
Pkg.activate(@__DIR__)

using ANTsRegistration
using Images, TestImages, Rotations, CoordinateTransformations, MappedArrays, FileIO, NRRD, Glob, Unitful, Statistics

img = testimage("cameraman")
tfm = Translation(125,250) ∘ LinearMap(RotMatrix(pi/50)) ∘ Translation(-125,-250)
img_rotated = warp(img, tfm)
fixed = img[50:300, 50:300]
moving = img_rotated[50:300, 50:300]

mreg, tforms = register(fixed, moving, Stage(fixed, Global("Rigid"), MI()))
@show tforms
