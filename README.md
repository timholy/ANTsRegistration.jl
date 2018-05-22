# ANTSRegistration

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

If you are registering an image sequence stored in NRRD format, and you see an error like

```sh
Description: itk::ERROR: NrrdImageIO(0x2b62880): ReadImageInformation: nrrd's #independent axes (3) doesn't match dimension of space in which orientation is defined (2); not currently handled
```

the remedy appears to be to delete the `space dimension` and `space
origin` fields from the header file. This is easier if you are using a
detached header file (`.nhdr`).

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

`stage` specifies that the metric (the way in which the two images are
compared) is `MI`, which stands for [mutual information](). (You can
specify various options such as the number of bins.) Alternatively you
can use `MeanSquares` or `CC`, which stands for neighborhood cross
correlation.

Finally, the final arguments in the example above signal that we want
to use a 4-level registration. For the first (coarsest) level, the
image will be shrunk by a factor of 8, smoothed over a 3-pixel radius,
and up to 1000 iterations will be used in trying to minimize the
metric. The main advantage of shrinking is performance, as small image
pairs require fewer pixelwise comparisons than large images.  The
advantage of smoothing is that you may increase the size of the
"attraction basin," making it easier to fall into a good minimum.

Once the rigid transformation has been found for this coarsest level,
it will be used to initialize the transformation for the next
level. The final level uses a `shrink` of 1 (meaning to use the images
at their provided size), a `smooth` of 0 (meaning no smoothing), and 5
iterations. This will ensure that the transformation doesn't miss
opportunities for sub-pixel alignment at the finest scale.


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

This will increase the odds that much of the registration will be an
affine alignment, and use warping only when necessary for higher
fidelity.


To register an image series in the file `movingfilename` to a single image `fixed`, use
```julia
motioncorr((infofilename, warpedfilename), fixed, movingfilename, pipeline)
```

See the help for more detail.
