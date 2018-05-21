module ANTSRegistration

using Images

export register, warp, Global, SyN, MeanSquares, MI, Stage

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
Global(mode::AbstractString) = Global(mode, 0.1)

Base.show(io::IO, st::Global) = print(io, st.mode, '[', st.gradientstep, ']')

struct SyN <: AbstractTransformation
    gradientstep::Float64
    updateFV::Float64
    totalFV::Float64
end
SyN() = SyN(0.1)
SyN(gradientstep) = SyN(gradientstep, 3, 0)

Base.show(io::IO, syn::SyN) = print(io, "SyN\[", syn.gradientstep, ',', syn.updateFV, ',', syn.totalFV, '\]')


abstract type AbstractMetric end

struct MeanSquares <: AbstractMetric
    metricweight::Float64
end
MeanSquares() = MeanSquares(1)
metric(cmd::Cmd, fixedname, movingname, m::MeanSquares) =
    `$cmd --metric MeanSquares\[$fixedname,$movingname,$(m.metricweight)\]`

struct MI <: AbstractMetric
    metricweight::Float64
    nhist::Int
    sampling::String
    percentage::Float64
end
MI() = MI(1, 32, "Regular", 0.25)
metric(cmd::Cmd, fixedname, movingname, m::MI) =
    `$cmd --metric MI\[$fixedname,$movingname,$(m.metricweight),$(m.nhist),$(m.sampling),$(m.percentage)\]`


function metric(io::IO, fixedname::AbstractString, movingname::AbstractString, m::AbstractMetric)
    print(io, metricstring(m), '[', fixedname, ',', movingname, ',')
    dumpparams(io, m)
    print(io, ']')
end
function metric(fixedname::AbstractString, movingname::AbstractString, m::AbstractMetric)
    io = IOBuffer()
    metric(io, fixedname, movingname, m)
    return String(take!(io))
end


struct Stage{N}
    transform::AbstractTransformation
    metric::AbstractMetric
    shrink::NTuple{N,Int}
    smoothing::NTuple{N,Int}
    iterations::NTuple{N,Int}
    convergencethresh::Float64
    convergencewindow::Int
end

Stage(img::AbstractArray, transform::AbstractTransformation, metric::AbstractMetric,
      shrink::NTuple{N,Int}, smooth::NTuple{N,Int}, iterations::NTuple{N,Int}) where N =
          Stage(transform, metric, shrink, smooth, iterations, default_convergence(img, transform)[2:3]...)
Stage(img::AbstractArray, transform::AbstractTransformation, metric::AbstractMetric,
      shrink::NTuple{N,Int}, smooth::NTuple{N,Int}) where N =
          Stage(transform, metric, shrink, smooth, default_convergence(img, transform)...)
Stage(img::AbstractArray, transform::AbstractTransformation, metric::AbstractMetric,
      shrink::NTuple{N,Int}) where N =
          Stage(img, transform, metric, shrink, default_smooth(img, transform))
Stage(img::AbstractArray, transform::AbstractTransformation, metric::AbstractMetric) =
    Stage(img, transform, metric, default_shrink(img, transform))
Stage(img::AbstractArray, transform::AbstractTransformation) =
    Stage(img, transform, MI())

function shcmd(cmd::Cmd, stage::Stage, fixedname, movingname)
    cmd = `$cmd --transform $(stage.transform)`
    cmd = metric(cmd, fixedname, movingname, stage.metric)
    cmd = `$cmd --convergence \[$(xstring(stage.iterations)),$(stage.convergencethresh),$(stage.convergencewindow)\]`
    cmd = `$cmd --shrink-factors $(xstring(stage.shrink)) --smoothing-sigmas $(xstring(stage.smoothing))vox`
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


function register(output, fixed::AbstractArray, moving::AbstractArray, pipeline::AbstractVector{<:Stage})
    maskfile = ""
    if any(isnan, fixed)
        # Create a mask
        error("not sure how a mask file should be implemented")
    end
    fixedname = write_nrrd(fixed)
    movingname = write_nrrd(moving)
    cmd = `$(ENV["ANTSPATH"])/antsRegistration -d $(sdims(fixed))`
    for pipe in pipeline
        cmd = shcmd(cmd, pipe, fixedname, movingname)
    end
    cmd = `$cmd --output `
    cmd = shcmd(cmd, output)
    @show cmd
    run(cmd)
end
function register(fixed::AbstractArray, moving::AbstractArray, pipeline::AbstractVector{<:Stage})
    up = userpath()
    outname = randstring(10)
    warpedname = joinpath(up, outname*".nrrd")
    output = (joinpath(up, outname*"_warp"), warpedname)
    register(output, fixed, moving, pipeline)
    return load(warpedname)
end

register(output, fixed::AbstractArray, moving::AbstractArray, pipeline::Stage) =
    register(output, fixed, moving, [pipeline])
register(fixed::AbstractArray, moving::AbstractArray, pipeline::Stage) =
    register(fixed, moving, [pipeline])


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

end # module
