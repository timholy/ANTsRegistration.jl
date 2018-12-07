module ANTsRegistration

using Images, Glob, Random, Unitful

export register, motioncorr, warp, Global, SyN, MeanSquares, CC, MI, Stage

# Per-user working path
function userpath()
    td = tempdir()
    user = ENV["USER"]
    up = joinpath(td, user, "ANTs")
    if !ispath(up)
        mkpath(up)
    end
    return up
end


abstract type AbstractTransformation end

const globaltransformations = ("Translation", "Rigid", "Similarity", "Affine")

struct Global <: AbstractTransformation
    mode::String
    gradientstep::Float64

    function Global(mode::AbstractString, gradientstep::Real)
        cmode = uppercase(mode[1])*lowercase(mode[2:end])
        cmode ∈ globaltransformations || error("mode $mode not recognized")
        return new(cmode, gradientstep)
    end
end
"""
    Global(mode, gradientstep=0.1)

Sets the mode of a global (affine) transformation. `mode` is one of the following strings:
  - "Translation"  (for image shifts only)
  - "Rigid"        (for shifts + rotations)
  - "Similarity"   (for rigid + scaling)
  - "Affine"       (for similarity + shear)
"""
Global(mode::AbstractString) = Global(mode, 0.1)

Base.show(io::IO, st::Global) = print(io, st.mode, '[', st.gradientstep, ']')

struct SyN <: AbstractTransformation
    gradientstep::Float64
    updateFV::Float64
    totalFV::Float64
end
"""
    Syn(gradientstep=0.1, updateFieldVariance=3, totalFieldVariance=0)

Sets the mode to diffeomorphic (warping) registration.
"""
SyN() = SyN(0.1)
SyN(gradientstep) = SyN(gradientstep, 3, 0)

Base.show(io::IO, syn::SyN) = print(io, "SyN[", syn.gradientstep, ',', syn.updateFV, ',', syn.totalFV, ']')


abstract type AbstractMetric end

struct MeanSquares <: AbstractMetric
    metricweight::Float64
end
"""
    MeanSquares(metricweight=1.0)

Compare images using mean-squared difference of voxel values.
"""
MeanSquares() = MeanSquares(1)
metric(cmd::Cmd, fixedname, movingname, m::MeanSquares) =
    `$cmd --metric MeanSquares\[$fixedname,$movingname,$(m.metricweight)\]`

struct MI <: AbstractMetric
    metricweight::Float64
    nhist::Int
    sampling::String
    percentage::Float64
end

"""
    MI(metricweight=1.0, nhist=32, sampling="Regular", percentage=0.25)

Compare images using mutual information.
"""
MI() = MI(1, 32, "Regular", 0.25)
metric(cmd::Cmd, fixedname, movingname, m::MI) =
    `$cmd --metric MI\[$fixedname,$movingname,$(m.metricweight),$(m.nhist),$(m.sampling),$(m.percentage)\]`

struct CC <: AbstractMetric
    metricweight::Float64
    radius::Int
end
"""
    CC(metricweight=1.0, radius=4)

Compare images using neighborhood cross correlation.
"""
CC() = CC(1, 4)
metric(cmd::Cmd, fixedname, movingname, m::CC) =
    `$cmd --metric CC\[$fixedname,$movingname,$(m.metricweight),$(m.radius)\]`


struct Stage{N,T}
    transform::AbstractTransformation
    metric::AbstractMetric
    shrink::NTuple{N,Int}
    smoothing::NTuple{N,T}
    iterations::NTuple{N,Int}
    convergencethresh::Float64
    convergencewindow::Int
end

Stage(img::AbstractArray, transform::AbstractTransformation, metric::AbstractMetric,
      shrink::NTuple{N,Int}, smooth::NTuple{N,Any}, iterations::NTuple{N,Int}) where N =
          Stage(transform, metric, shrink, smooth, iterations, default_convergence(img, transform)[2:3]...)
Stage(img::AbstractArray, transform::AbstractTransformation, metric::AbstractMetric,
      shrink::NTuple{N,Int}, smooth::NTuple{N,Any}) where N =
          Stage(transform, metric, shrink, smooth, default_convergence(img, transform)...)
Stage(img::AbstractArray, transform::AbstractTransformation, metric::AbstractMetric,
      shrink::NTuple{N,Int}) where N =
          Stage(img, transform, metric, shrink, default_smooth(img, transform))
Stage(img::AbstractArray, transform::AbstractTransformation, metric::AbstractMetric) =
    Stage(img, transform, metric, default_shrink(img, transform))

"""
    stage = Stage(fixedimg, transform, metric=MI(), shrink=(8,4,2,1), smooth=(3,2,1,0), iterations=(1000,500,250,0))

Create a registration stage. `transform` is an AbstractTransformation;
currently [`Global`](@ref) and [`Syn`](@ref) are supported. `metric`
is one of [`MI`](@ref), [`MeanSquares`](@ref), or [`CC`](@ref). The
remaining arguments are tuple arguments for multilevel optimization;
each corresponding entry indicates the shrinking factor (8 means 8×
reduced in size), smoothing size (in voxels), and number of iterations
allowed when optimizing.
"""
Stage(img::AbstractArray, transform::AbstractTransformation) =
    Stage(img, transform, MI())

function shcmd(cmd::Cmd, stage::Stage, fixedname, movingname; ismc::Bool=false)
    cmd = `$cmd --transform $(stage.transform)`
    cmd = metric(cmd, fixedname, movingname, stage.metric)
    if ismc
        cmd = `$cmd --iterations $(xstring(stage.iterations))`
        cmd = `$cmd --shrinkFactors $(xstring(stage.shrink)) --smoothingSigmas $(xstring(stage.smoothing))`
    else
        cmd = `$cmd --convergence \[$(xstring(stage.iterations)),$(stage.convergencethresh),$(stage.convergencewindow)\]`
        cmd = `$cmd --shrink-factors $(xstring(stage.shrink)) --smoothing-sigmas $(xqstring(stage.smoothing))$(ustring(stage.smoothing))`
    end
    return cmd
end

default_shrink(img::AbstractArray, transform) = default_shrink(size_spatial(img), transform)
function default_shrink(sz::Dims, transform::Global)
    islarge = any(s->s>=256, sz)
    return islarge ? (12,8,4,2) : (8,4,2,1)
end

function default_shrink(sz::Dims, transform::SyN)
    islarge = any(s->s>=256, sz)
    return islarge ? (10,6,4,2,1) : (8,4,2,1)
end


default_smooth(img::AbstractArray, transform) = default_smooth(size_spatial(img), transform)
function default_smooth(sz::Dims, transform::Global)
    islarge = any(s->s>=256, sz)
    return islarge ? (4,3,2,1) : (3,2,1,0)
end

function default_smooth(sz::Dims, transform::SyN)
    islarge = any(s->s>=256, sz)
    return islarge ? (5,3,2,1,0) : (3,2,1,0)
end


default_convergence(img::AbstractArray, transform) = default_convergence(size_spatial(img), transform)
function default_convergence(sz::Dims, transform::Global)
    islarge = any(s->s>=256, sz)
    return ((1000,500,250,0), 1e-6, 10)
end

function default_convergence(sz::Dims, transform::SyN)
    islarge = any(s->s>=256, sz)
    return islarge ? ((100,100,70,50,0), 1e-6, 10) : ((100,70,50,0), 1e-6, 10)
end


function register(output, nd::Int, fixedname::AbstractString, movingname::AbstractString, pipeline::AbstractVector{<:Stage}; histmatch::Bool=false, winsorize=nothing, verbose::Bool=false)
    cmd = `$(ENV["ANTSPATH"])/antsRegistration -d $nd`
    if verbose
        cmd = `$cmd -v`
    end
    if histmatch
        cmd = `$cmd --use-histogram-matching`
    end
    if isa(winsorize, Tuple{Real,Real})
        cmd = `$cmd --winsorize-image-intensities \[$(winsorize[1]),$(winsorize[2])\]`
    end
    for pipe in pipeline
        cmd = shcmd(cmd, pipe, fixedname, movingname)
    end
    cmd = `$cmd --output `
    cmd = shcmd(cmd, output)
    if verbose
        @show cmd
    end
    run(cmd)
end

function register(output, fixed::AbstractArray, moving::AbstractArray, pipeline::AbstractVector{<:Stage}; kwargs...)
    maskfile = ""
    if any(isnan, fixed)
        # Create a mask
        error("not sure how a mask file should be implemented")
    end
    fixedname = write_nrrd(fixed)
    movingname = write_nrrd(moving)
    imgw = register(output, sdims(fixed), fixedname, movingname, pipeline; kwargs...)
    rm(movingname)
    rm(fixedname)
end

function register(fixed::AbstractArray, moving, pipeline::AbstractVector{<:Stage}; kwargs...)
    up = userpath()
    outname = randstring(10)
    tfmname = joinpath(up, outname*"_warp")
    warpedname = joinpath(up, outname*".nrrd")
    output = (tfmname, warpedname)
    register(output, fixed, moving, pipeline; kwargs...)
    imgw = load(warpedname)
    rm(warpedname)
    for tfmfile in glob(outname*"_warp"*"*.mat", up)
        rm(tfmfile)
    end
    for tfmfile in glob(outname*"_warp"*"*.nii.gz", up)
        rm(tfmfile)
    end
    return imgw
end

register(output, fixed::AbstractArray, moving, pipeline::Stage; kwargs...) =
    register(output, fixed, moving, [pipeline]; kwargs...)

"""
    imgw = register(fixed, moving, pipeline; kwargs...)

Return a version of `moving` that has been warped to match
`fixed`. `pipeline` is a single [`Stage`](@ref) or a vector of stages.

# Example

    ## Create some images
    using Images, TestImages, Rotations, CoordinateTransformations
    img = testimage("cameraman")
    tfm = Translation(125,250) ∘ LinearMap(RotMatrix(pi/50)) ∘ Translation(-125,-250)
    img_rotated = warp(img, tfm)
    fixed = img[50:300, 50:300]
    moving = img_rotated[50:300, 50:300]

    ## Perform the registration
    rigid = Stage(fixed, Global("Rigid"))
    syn = Stage(fixed, SyN())
    imgw = register(fixed, moving, [rigid,syn])
"""
register(fixed::AbstractArray, moving, pipeline::Stage; kwargs...) =
    register(fixed, moving, [pipeline]; kwargs...)


"""
    motioncorr(output, fixed, movingname, pipeline; kwargs...)

Perform motion correction in an image series. All images are
registered to `fixed`.  `movingname` is the name of the file that
contains all the individual images, and must be stored in a format
that ITK can handle (see the README). See [`register`](@ref) for more
information about the other arguments.

`output` is generally a tuple of the form `(deformationinfoprefix,
warpedimagefilename)`.

# Example

    output = (joinpath(mydir, "seriesinfo"), joinpath(mydir, "series_warped.nrrd"))
    # Set the number of iterations to be pretty modest (this example uses 2 levels)
    stage = Stage(fixed, Global("Rigid"), CC(), (2,1), (1,0), (15,3))
    motioncorr(output, fixed, movingfilename, stage; verbose=true)
"""
function motioncorr(output, fixed::AbstractArray, movingname::AbstractString, pipeline::AbstractVector{<:Stage}; kwargs...)
    maskfile = ""
    if any(isnan, fixed)
        # Create a mask
        error("not sure how a mask file should be implemented")
    end
    fixedname = write_nrrd(fixed)
    motioncorr(output, sdims(fixed), fixedname, movingname, pipeline; kwargs...)
    rm(fixedname)
end
motioncorr(output, fixed::AbstractArray, movingname::AbstractString, stage::Stage; kwargs...) =
    motioncorr(output, fixed, movingname, [stage]; kwargs...)

function motioncorr(output, nd::Int, fixedname::AbstractString, movingname::AbstractString, pipeline::AbstractVector{<:Stage}; verbose::Bool=false)
    cmd = `$(ENV["ANTSPATH"])/antsMotionCorr -u 1 -e 1 -d $nd`
    if verbose
        cmd = `$cmd -v`
    end
    for pipe in pipeline
        cmd = shcmd(cmd, pipe, fixedname, movingname; ismc=true)
    end
    cmd = `$cmd --output `
    cmd = shcmd(cmd, output)
    if verbose
        @show cmd
    end
    run(cmd)
end


shcmd(cmd::Cmd, str::AbstractString) = `$cmd $str`
function shcmd(cmd::Cmd, strs::Tuple{Vararg{S}}) where S <: AbstractString
    io = IOBuffer()
    print(io, strs[1])
    for i = 2:length(strs)
        print(io, ',', strs[i])
    end
    cmdstr = String(take!(io))
    return `$cmd \[$cmdstr\]`
end

function write_nrrd(img::AbstractArray)
    fn = joinpath(userpath(), randstring(10)*".nrrd")
    save(fn, img)
    return fn
end

function xstring(io::IO, t::Tuple)
    print(io, t[1])
    for i = 2:length(t)
        print(io, 'x', t[i])
    end
end

function xstring(t::Tuple)
    io = IOBuffer()
    xstring(io, t)
    return String(take!(io))
end

xqstring(t) = xstring(to_mm_or_vox.(t))

to_mm_or_vox(x::Number) = x
# to_mm_or_vox(x::Quantity) = float(ustrip(uconvert(u"mm", x)))
to_mm_or_vox(x::Quantity) = float(ustrip(x))  # it seems that ITK treats any physical units as "mm"

ustring(t) = ustring(t[1])
ustring(x::Number) = "vox"
ustring(x::Quantity) = "mm"   # seems to require mm

end # module
