using ANTsRegistration, Images, TestImages, Rotations, CoordinateTransformations, MappedArrays, FileIO, NRRD, Glob
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

@testset "Options" begin
    fixed = [0 0 0 0 0 0 0;
             0 0 0 0 0 0 0;
             0 0 10 0 0 0 0;
             0 0 0 0 0 0 0;
             0 0 0 0 1 1 0;
             0 0 0 0 0 0 0]
    moving = [11 0 0 0 0 0 0;
              0 0 0 0 0 0 0;
              0 0 8 8 0 0 0;
              zeros(3, 7)]
    stage = Stage(fixed, Global("Translation"), MeanSquares(), (1,1), (2,0), (1000,1000))
    imgw = register(fixed, moving, stage)
    @test_broken sum(imgw) > 0
    @test_broken sum(mappedarray(diff2_0, fixed, histmatch(imgw, fixed, 12))) > 1e-3*sum(abs2f, fixed)
    imgw = register(fixed, moving, stage; histmatch=true)
    @test_broken sum(imgw) > 0
    @test sum(mappedarray(diff2_0, fixed, histmatch(imgw, fixed, 12))) < 1e-3*sum(abs2f, fixed)
    fixed = rand(100,100)
    moving = copy(fixed)
    fixed[50,50] = 1000
    moving[50,51] = 1000
    imgw = register(fixed, moving, stage)
    @test median(abs.(imgw[:,1:end-1] - moving[:,2:end])) < 1e-3
    imgw = register(fixed, moving, stage; winsorize=(0.05,0.95))
    @test median(abs.(imgw - moving)) < 0.01
end

@testset "TimeSeries" begin
    img = testimage("cameraman")
    fixed = img[50:300, 50:300]
    ## Create a movie with a time dimension and progressive rigid transformation
    testdir = ANTsRegistration.userpath()
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
    motioncorr(output, fixed, hdrfile, Stage(fixed, Global("Rigid"), CC(), (2,1), (1,0), (15,3)))
    imgw = load(output[2])
    for i = 1:size(imgw, 3)
        @test sum(mappedarray(diff2_0, fixed, imgw[:,:,i]/255)) < 1e-3*sum(abs2f, fixed)
    end
    for tfmfile in glob("seriesinfo*.csv", testdir)
        rm(tfmfile)
    end
    rm(output[2])

    # Also do it via register. Exercise the low-level call that uses only filenames.
    moving = load(hdrfile)
    fixedname = joinpath(testdir, "series_fixed.png")
    save(fixedname, fixed)
    movingtmp = joinpath(testdir, "series_frame.png")
    stageMS = Stage(fixed, Global("Rigid"), MeanSquares())
    stageCC = Stage(fixed, Global("Rigid"), CC())
    for i = 2:size(moving, 3)
        save(movingtmp, moving[:,:,i])
        output = (joinpath(testdir, "info$i"), joinpath(testdir, "series_warp.nrrd"))
        register(output, sdims(fixed), fixedname, movingtmp, [stageMS])
        imgw = load(output[2])'/255  # FIXME NRRD orientation
        @test sum(mappedarray(diff2_0, fixed, imgw)) < 1e-3*sum(abs2f∘Float64, fixed)
        # register(output, sdims(fixed), fixedname, movingtmp, [stageCC])
        # imgw = load(output[2])/255
        # @test sum(mappedarray(diff2_0, fixed, imgw)) < 1e-3*sum(abs2f, fixed)
    end
    for tfmfile in glob("info*.mat", testdir)
        rm(tfmfile)
    end
    rm(output[2])
    rm(rawfile)
    rm(hdrfile)
    rm(fixedname)
    rm(movingtmp)
end
