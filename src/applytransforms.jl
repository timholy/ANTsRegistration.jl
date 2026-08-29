abstract type AbstractAntsInterpolation end

####Interpolation mode
interpolation_mode = ("Linear", "NearestNeighbor", "MultiLabel", "Gaussian", "BSpline", "CosineWindowedSinc", "HammingWindowedSinc", "LanczosWindowedSinc", "GenericLabel") 

struct Linear<:AbstractAntsInterpolation
    mode::AbstractString
end
Linear() = Linear("Linear")

struct NearestNeighbor<:AbstractAntsInterpolation
    mode::AbstractString
end
NearestNeighbor() = NearestNeighbor("NearestNeighbor")

### FIXME multilabel
struct MultiLabel{N}<:AbstractAntsInterpolation
   mode::AbstractString
   sigma::NTuple{N, Float64}
   alpha::Float64
end

### FIXME gaussian
struct Gaussian{N}<:AbstractAntsInterpolation
   mode::AbstractString
   sigma::NTuple{N, Number}
   alpha::Float64
end

struct BSpline<:AbstractAntsInterpolation
    mode::AbstractString
    order::Int
end
BSpline() = BSpline("BSpline", 3)
BSpline(order) = BSpline("BSpline", order)

struct CosineWindowedSinc<:AbstractAntsInterpolation
    mode::AbstractString
end
CosineWindowedSinc() = CosineWindowedSinc("CosineWindowedSinc")

struct WelchWindowedSinc<:AbstractAntsInterpolation
    mode::AbstractString
end
WelchWindowedSinc() = WelchWindowedSinc("WelchWindowedSinc")

struct HammingWindowedSinc<:AbstractAntsInterpolation
    mode::AbstractString
end
HammingWindowedSinc() = HammingWindowedSinc("HammingWindowedSinc")

struct LanczosWindowedSinc<:AbstractAntsInterpolation
    mode::AbstractString
end
LanczosWindowedSinc() = LanczosWindowedSinc("LanczosWindowedSinc")

struct GenericLabel<:AbstractAntsInterpolation
    mode::AbstractString
    interpolator::AbstractString
end
GenericLabel() = GenericLabel("GenericLabel", "Linear")
GenericLabel(interpolator) = GenericLabel("GenericLabel", interpolator)

"""ITKTransform struct."""
struct ITKTransform
    mode::AbstractString
    version::AbstractString
    tag::AbstractString #Not sure what it means in ANTs.
    transform::AbstractString
    parameters::NTuple
    fixedparameters::NTuple
end

"""Transformation setup for antsApplyTransform"""
struct Tform
    transform::ITKTransform
    useInverse::Int
end
Tform(transform::ITKTransform) = Tform(transform, 0)
Tform(transformFileName::AbstractString) = Tform(load_itktform(transformFileName))
Tform(transformFileName::AbstractString, useInverse::Int) = Tform(load_itktform(transformFileName), useInverse)

"""Point.
`x` is the row position from `Images.imshow`.
`y` is the column position from `Images.imshow`.
If `z` and `t` are not used, set them as 0.
"""
struct Point
    x::Number
    y::Number
    z::Number
    t::Number
end
Point(x::Number, y::Number, z::Number) = Point(x, y, z, 0)
Point(x::Number, y::Number) = Point(x, y, 0, 0)

function get_tempname(tag::AbstractString)
    tfmname = joinpath(userpath(), randstring(10)*tag) #temporary transform file names
    return tfmname
end
get_tempname(noutputs::Int, tag::AbstractString) = [get_tempname(tag) for i in 1:noutputs]
get_tempname(noutputs::Int) = [get_tempname() for i in 1:noutputs]
get_tempname() = get_tempname("")

function applyTransforms(outputFileName, nd::Int, tforms::Vector{Tform}, referenceFileName::AbstractString, inputFileName::AbstractString; interpolation::AbstractAntsInterpolation = Linear(), input_imagetype = 0, output_datatype = "default", float::Bool = false, default_value = missing, verbose::Bool=false, suppressout::Bool=true)
    tfmnames = get_tempname(length(tforms), "_tfm.txt");
    cmd = `antsApplyTransforms -o $outputFileName -d $nd -r $referenceFileName -i $inputFileName --input-image-type $input_imagetype --output-data-type $output_datatype`
    # Add interpolation method
    if any(x -> isa(interpolation, x), [Linear, NearestNeighbor, CosineWindowedSinc, WelchWindowedSinc, HammingWindowedSinc, LanczosWindowedSinc])  
        cmd = `$cmd -n $(interpolation.mode)`
    elseif any(x -> isa(interpolation, x), [MultiLabel, Gaussian])
        cmd = `$cmd -n \[$(interpolation.mode), $(interpolation.sigma), $(interpolation.alpha)\]`
    elseif any(x -> isa(interpolation, x), [BSpline])
        cmd = `$cmd -n \[$(interpolation.mode), $(interpolation.order)\]`
    end
    # Add transformations
    for (i, tform) in enumerate(tforms)
        save_itktform(tfmnames[i], tform.transform)
        cmd = `$cmd -t \[$(tfmnames[i]), $(tform.useInverse)\]`
    end
    # Other options
    if verbose
        cmd = `$cmd -v 1`
    end
    if !ismissing(default_value)
        cmd = `$cmd --default-value $default_value`
    end
    # output display
    if verbose
        @show cmd
    end
    if suppressout
        @suppress_out run(cmd)
    else
        run(cmd)
    end
end

function applyTransforms(outputFileName, tforms::Vector{Tform}, reference::AbstractArray, input::AbstractArray; kwargs...)
    refimg = write_nrrd(reference);
    inputimg = write_nrrd(input);
    applyTransforms(outputFileName, sdims(reference), tforms, refimg, inputimg; kwargs...)
    rm(refimg)
    rm(inputimg)
end

""" Apply Transforms to an image and returns a transformed image.

`tforms` is a vector of `Tform`. This requires ITKTransform, most likely obtained from the `register` function.
A `reference` image defines spacing, origin, size, and the direction of the output warped image.
Transformation is applied to an `input` image.

See antsApplyTransforms.
"""
function applyTransforms(tforms::Vector{Tform}, reference::AbstractArray, input::AbstractArray; kwargs...)
    outputFileName = get_tempname("_warp.nrrd")
    refimg = write_nrrd(reference);
    inputimg = write_nrrd(input);
    applyTransforms(outputFileName, sdims(reference), tforms, refimg, inputimg; kwargs...)
    imgw = load(outputFileName)
    rm(refimg)
    rm(inputimg)
    rm(outputFileName)
    return imgw
end

function applyTransformsToPoints(outputFileName::AbstractString, nd::Int, tforms::Vector{Tform}, inputFileName::AbstractString; precision::Bool = false)
    tfmnames = get_tempname(length(tforms), "_tfm.txt")
    cmd = `antsApplyTransformsToPoints -o $outputFileName -d $nd -i $inputFileName`
    if precision
        cmd = `$cmd --precision 1`
    end
    for (i, tform) in enumerate(tforms)
        save_itktform(tfmnames[i], tform.transform)
        cmd = `$cmd -t \[$(tfmnames[i]), $(tform.useInverse)\]`
    end
    run(cmd)
end

function applyTransformsToPoints(outputFileName::AbstractString, nd::Int, tforms::Vector{Tform}, points::DataFrame; precision::Bool = false)
    tmpinputname = get_tempname("_point.csv")
    CSV.write(tmpinputname, points)
    applyTransformsToPoints(outputFileName, nd, tforms, tmpinputname; precision = precision)
    df_tformed = CSV.read(outputFileName, DataFrame)
    rm(tmpinputname)
    return df_tformed
end

function applyTransformsToPoints(outputFileName::AbstractString, nd::Int, tforms::Vector{Tform}, points::Vector{Point}; precision::Bool = false)
    x = map(p -> p.x, points)
    y = map(p -> p.y, points)
    z = map(p -> p.z, points)
    t = map(p -> p.t, points)
    df = DataFrame(x = x, y = y, z = z, t = t) #Make data frames
    df_tformed = applyTransformsToPoints(outputFileName, nd, tforms, df; precision = precision)
    points_tformed = map(p -> Point(p...), zip(df_tformed.x, df_tformed.y, df_tformed.z, df_tformed.t))
    return points_tformed 
end

function applyTransformsToPoints(nd::Int, tforms::Vector{Tform}, points::DataFrame; precision::Bool = false)
    tmpoutname = get_tempname("_point.csv")
    df_tformed = applyTransformsToPoints(tmpoutname, nd, tforms, points; precision = precision)
    rm(tmpoutname)
    return df_tformed
end

"""Apply Transforms to point
`nd` is dimensionality (2/3/4).
`tforms` is a vector of `Tform`. `Tform` requires `ITKTransform`, most likely obtained from the `register` function.
`points` are a vector of `Point` with x, y, z, and t. Instead of `Point`, a data frame also works.

See `antsApplyTransformsToPoint` from ANTs.
"""
function applyTransformsToPoints(nd::Int, tforms::Vector{Tform}, points::Vector{Point}; precision::Bool = false)
    tmpoutname = get_tempname("_point.csv")
    points_tformed = applyTransformsToPoints(tmpoutname, nd, tforms, points; precision = precision)
    rm(tmpoutname)
    return points_tformed
end

