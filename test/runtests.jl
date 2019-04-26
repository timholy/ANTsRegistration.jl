using ANTsRegistration, Images, TestImages, Rotations, CoordinateTransformations, MappedArrays, FileIO, NRRD, Glob
using Unitful
using Test, Statistics

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
    # Histogram matching
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
    pfixed  = padarray(fixed, Fill(0, (3,3))).parent
    pmoving = padarray(moving, Fill(0, (3,3))).parent
    stage = Stage(pfixed, Global("Translation"), MeanSquares(), (1,1), (3,0), (1000,1000))
    imgw = register(pfixed, pmoving, stage)
    @test sum(mappedarray(diff2_0, pfixed, adjust_histogram(Matching(), imgw, pfixed, 256))) > 1e-3*sum(abs2f, fixed)
    imgw = register(pfixed, pmoving, stage; histmatch=true)
    @test sum(mappedarray(diff2_0, pfixed, adjust_histogram(Matching(), imgw, pfixed, 0:12))) < 1e-3*sum(abs2f, fixed)

    # Winsorizing
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
    header = NRRD.headerinfo(N0f8, (Axis{:y}(axes(fixed, 1)), Axis{:x}(axes(fixed, 2)), Axis{:time}(axes(θs, 1))))
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
    imgw = moving = nothing; GC.gc(); # Unlink variables from temp files. added to prevent IOerror when run on Windows. 
    rm(output[2])
    rm(rawfile)
    rm(hdrfile)
    rm(fixedname)
    rm(movingtmp)
end

@testset "Anisotropic spacing" begin
    img = testimage("cameraman")
    fixed0 = img[50:300, 50:300]
    tfm = Translation(250,250) ∘ LinearMap(RotMatrix(-pi/4)) ∘ Translation(-250,-250)
    img_rotated = warp(img, tfm)
    moving0 = img_rotated[75:320, 75:325]
    ps = (1u"mm", 2u"mm")
    # If you don't trim the edges, the mixing with zeros "anchors" the image in place.
    # see resolution of https://github.com/ANTsX/ANTs/issues/675
    fixed = AxisArray(restrict(fixed0, 2)[:,2:end-1], (:y, :x), ps)
    moving = AxisArray(restrict(moving0, 2)[:,2:end-1], (:y, :x), ps)
    stage = Stage(fixed, Global("Rigid"), MeanSquares(), (1,), (25u"mm",), (1000,))
    # Both of the following should register well: fixed0 and moving0 are isotropic, so
    # the rotation is a rotation; fixed and moving encode their pixel spacing, so even
    # though "squashed" horizontally they also differ by a rotation (once the pixel spacing is accounted for)
    imgw0 = register(fixed0, moving0, stage)/255
    imgw = AxisArray(register(fixed, moving, stage), (:y, :x), ps)
    @test 2*mean(mappedarray(diff2_0, fixed0, imgw0)) >= mean(mappedarray(diff2_0, fixed, imgw))
    # In contrast, if we take away the pixel spacing info, we should not register well
    imgws = register(fixed.data, moving.data, stage)
    @test mean(mappedarray(diff2_0, fixed.data, imgws)) > 3*mean(mappedarray(diff2_0, fixed, imgw))
end
