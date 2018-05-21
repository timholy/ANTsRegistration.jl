using ANTSRegistration, Images, TestImages, Rotations, CoordinateTransformations, MappedArrays
using Base.Test

diff2_0(f, m) = m == 0 ? 0.0 : (Float64(f) - Float64(m))^2
abs2f(x) = abs2(Float64(x))

@testset "Single stage" begin
    img = testimage("cameraman")
    tfm = Translation(125,250) ∘ LinearMap(RotMatrix(pi/50)) ∘ Translation(-125,-250)
    img_rotated = warp(img, tfm)
    fixed = img[50:300, 50:300]
    moving = img_rotated[50:300, 50:300]
    mreg = register(fixed, moving, Stage(fixed, Global("Rigid")))/255
    @test sum(mappedarray(diff2_0, fixed, mreg)) < 1e-3*sum(abs2f, fixed)

    # For a bigger rotation we need to customize
    tfm = Translation(125,250) ∘ LinearMap(RotMatrix(pi/12)) ∘ Translation(-125,-250)
    img_rotated = warp(img, tfm)
    fixed = img[50:300, 50:300]
    moving = img_rotated[50:300, 50:300]
    mreg = register(fixed, moving, Stage(fixed, Global("Rigid"), MeanSquares(), (8,4,2,1), (10,7,3,0)))/255
    @test sum(mappedarray(diff2_0, fixed, mreg)) < 1e-3*sum(abs2f, fixed)
end
