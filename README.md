# ANTsRegistration

[![Build Status](https://travis-ci.org/timholy/ANTSRegistration.jl.svg?branch=master)](https://travis-ci.org/timholy/ANTSRegistration.jl)

[![codecov.io](http://codecov.io/github/timholy/ANTSRegistration.jl/coverage.svg?branch=master)](http://codecov.io/github/timholy/ANTSRegistration.jl?branch=master)

This provides a Julia wrapper around the
[Advanced Normalization Tools](https://stnava.github.io/ANTs/) image
registration and motion correction suite.

## Installation
```julia
] add https://github.com/timholy/ANTsRegistration.jl
```

## Usage

### Image data and file format

If you are passing the data via filenames, ensure that you have stored
your images in an ITK-readable
format. [NRRD](https://github.com/JuliaIO/NRRD.jl) is recommended. For
those performing acquisition with Imagine, you can write out an
[NRRD header](https://github.com/timholy/ImagineFormat.jl#converting-to-nrrd).

Unfortunately, [the NRRD format](http://teem.sourceforge.net/nrrd/format.html)
lacks a file-validator, and a few aspects of the standard description seem
to leave room for interpretation. If you encounter bugs, it is possible that
Julia and ITK differ with respect to the implementation of the NRRD header.
Try copying the command and running it the shell prompt, then report the error
[here](https://github.com/JuliaIO/NRRD.jl/issues/new).

### Performing registration

#### Stages

The process of registering images can be broken down into stages, and multiple stages can be cascaded together. A stage is specified as follows:

```julia
stage = Stage(fixedimg, transform, metric=MI(), shrink=(8,4,2,1), smooth=(3,2,1,0), iterations=(1000,500,250,5))
```

The transform might be one of the following:
```julia
transform = Global("Rigid")
transform = Global("Affine")
transform = Syn()
```
The last one is for a diffeomorphism (warping) registration.

This particular `stage` uses the `MI` metric for comparing the two
images.  `MI` is short for
[mutual information](https://en.wikipedia.org/wiki/Mutual_information),
a generalization of the notion of cross-correlation.  This can be a
good choice particularly when the images differ in ways other than
just a spatial transformation, for example when they may be collected
by different imaging modalities or exhibit intensity differences due
to calcium transients. (With `MI` you can optionally specify various
parameters such as the number of histogram bins.) Alternatively you
can use `MeanSquares` (where the images are compared based on their
mean-squared-difference) or `CC` (which stands for neighborhood cross
correlation).

Finally, the last arguments in the example above indicate that we want
to use a 4-level registration. For the first (coarsest) level, the
image will be shrunk by a factor of 8, smoothed over a 3-pixel radius,
and then aligned, allowing the parameters to be tweaked up to 1000
times when trying to minimize the metric. Choosing to shrink can
improve performance, because small image pairs require fewer pixelwise
comparisons than large images, as long as you don't shrink so much
that features useful for alignment are eliminated.  Likewise,
smoothing can help find a good minimum by increasing the size of the
"attraction basin," as long as you don't blur out sharp features that
actually aid alignment.

Once the rigid transformation has been found for this coarsest level,
it will be used to initialize the transformation for the next
level. The final level uses a `shrink` of 1 (meaning to use the images
at their provided size), a `smooth` of 0 (meaning no smoothing), and 5
iterations. This will ensure that the transformation doesn't miss
opportunities for sub-pixel alignment at the finest scale.

All parameters after `transform` have default values, so you only need
to assign them if you need to control them more precisely.

**Note on physical units**: if your images have anisotropic resolution,
you should strongly consider using physical units for your smoothing.
For example,

```julia
using Unitful: μm
smooth=(50μm,5μm)
```

would be appropriate for a two-iteration stage.

#### Top-level API

To register the single image `moving` to the single image `fixed`, use
```julia
imgw = register(fixed, moving, pipeline; kwargs...)
```

where `pipeline` is a single `Stage` or a vector of stages. For
example, you can begin with an affine registration followed by a
deformable registration:

```julia
stageaff = Stage(fixed, Global("Affine"))
stagesyn = Stage(fixed, Syn())
imgw = register(fixed, moving, [stageaff,stagesyn]; kwargs...)
```

This choice will align the images as well as possible (given the
default parameters) using a pure-affine transformation, and then
introduce warping where needed to improve the alignment.

If instead you'd like to correct for motion in an image sequence, consider
```julia
motioncorr((infofilename, warpedfilename), fixed, movingfilename, pipeline)
```
Here you represent the moving image sequence via its filename. The
first argument stores the names of the files to which the data should
be written. Of course, you can alternatively call `register` iteratively
for each image in the series.

For more detailed information, see the help on individual types and
functions.

### Some notes for working on Windows
Currently, ANTs does not actively support Windows system, officially, they suggest using Linux subsystem instead.

And [the latest binaries they built and released](https://github.com/ANTsX/ANTs/releases/tag/v2.1.0) is 4 years old, and it may not work on every machine. Because of this, if you just install this package on Windows system and use the default Binary files downloaded by `BinaryProvider` in `deps.jl`, [it's highly possible it will not work](https://github.com/ANTsX/ANTs/issues/339).

However, it's tested ANTs can be built from source in Windows system. (Windows 10, VS 2017) (Although officially they suggest [compiling under Linux subsystem](https://github.com/ANTsX/ANTs/wiki/Compiling-ANTs-on-Windows-10) ). Conceptually it's the same as compiling under MacOS / Linux. It majorly involves using `Cmake` to determine the specificities of your machine / system and automatically generate the configuration file. And then using the compilers of Windows system (e.g. Visual Studio) to build the binary files from source.

After building from source (usually takes hours), you get your path containing the binary files (supposing it's `D:\ANTs_2.1.0_Windows_new_build\bin\Release`). Then, we have to inform `ANTsRegistration.jl` package of the correct location of the binary files. Currently, we have to do it by manually changing `deps.jl` file generated when building this julia package.

For example, change the relevent lines in `deps.jl` to be
```julia
const ants_bin_dir = raw"D:\ANTs_2.1.0_Windows_new_build\bin\Release\\" # use this variable for outside ANTs binary
const ants = joinpath(ants_bin_dir, "ANTS.exe") #joinpath(dirname(@__FILE__), "usr\\bin\\ANTS.exe")
const antsRegistration = joinpath(ants_bin_dir, "antsRegistration.exe") # joinpath(dirname(@__FILE__), "usr\\bin\\antsRegistration.exe")
const antsMotionCorr = joinpath(ants_bin_dir, "antsMotionCorr.exe")
```
Then, hopefully, the package can pass all of the tests!
