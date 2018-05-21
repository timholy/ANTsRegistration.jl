using ANTSRegistration, Images, TestImages, Rotations, CoordinateTransformations, MappedArrays, FileIO, NRRD
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

@testset "Multistage" begin
    img = testimage("cameraman")
    tfm = Translation(125,250) ∘ LinearMap(RotMatrix(pi/50)) ∘ Translation(-125,-250)
    img_rotated = warp(img, tfm)
    fixed = img[50:300, 50:300]
    moving = img_rotated[50:300, 50:300]
    rigid = Stage(fixed, Global("Rigid"))
    syn = Stage(fixed, SyN())
    mreg = register(fixed, moving, [rigid,syn])/255
    @test sum(mappedarray(diff2_0, fixed, mreg)) < 1e-3*sum(abs2f, fixed)
end

@testset "TimeSeries" begin
    img = testimage("cameraman")
    fixed = img[50:300, 50:300]
    ## Create a movie with a time dimension and progressive rigid transformation
    testdir = ANTSRegistration.userpath()
    rawfile = joinpath(testdir, "series.raw")
    θs = ((0:5)/5)*(π/50)
    open(rawfile, "w") do io
        for θ in θs
            tfm = Translation(125,250) ∘ LinearMap(RotMatrix(θ)) ∘ Translation(-125,-250)
            img_rotated = warp(img, tfm)
            moving = img_rotated[50:300, 50:300]
            write(io, moving)
        end
    end
    hdrfile = joinpath(testdir, "series.nhdr")
    header = NRRD.headerinfo(N0f8, (Axis{:y}(indices(fixed, 1)), Axis{:x}(indices(fixed, 2)), Axis{:time}(indices(θs, 1))))
    header["datafile"] = rawfile
    # Hack to overcome ITK's limitations in reading NRRD files
    delete!(header, "space dimension")
    delete!(header, "space origin")
    open(hdrfile, "w") do io
        write(io, magic(format"NRRD"))
        NRRD.write_header(io, "0004", header)
    end

    ## Register the sequence
    output = (joinpath(testdir, "seriesinfo"), joinpath(testdir, "series_warped.nrrd"))
    motioncorr(output, fixed, hdrfile, Stage(fixed, Global("Rigid"), CC(), (2,1), (1,0), (15,3)); verbose=true)
    imgw = load(output[2])
    for i = 1:size(imgw, 3)
        @test sum(mappedarray(diff2_0, fixed, imgw[:,:,i]/255)) < 1e-3*sum(abs2f, fixed)
    end
end
