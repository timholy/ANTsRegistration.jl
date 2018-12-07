# ANTsRegistration

[![Build Status](https://travis-ci.org/timholy/ANTSRegistration.jl.svg?branch=master)](https://travis-ci.org/timholy/ANTSRegistration.jl)

[![codecov.io](http://codecov.io/github/timholy/ANTSRegistration.jl/coverage.svg?branch=master)](http://codecov.io/github/timholy/ANTSRegistration.jl?branch=master)

This provides a Julia wrapper around the
[Advanced Normalization Tools](https://stnava.github.io/ANTs/) image
registration and motion correction suite.

## Installation

To use this package you need to install ANTS manually and define the
`ANTSPATH` variable. For example, I built
ANTS
[from source](https://brianavants.wordpress.com/2012/04/13/updated-ants-compile-instructions-april-12-2012/)
and then added

```sh
export ANTSPATH="/home/tim/src/antsbin/bin"
```

to my `.bashrc` file. If this code throws a `Key Error`, then most
likely you didn't define this variable or need to execute `source
~/.bashrc` before launching Julia.

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
