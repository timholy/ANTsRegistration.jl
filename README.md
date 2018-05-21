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

First, if you are correcting for image motion in a timeseries, make
sure your input image has been properly annotated as an AxisArray with
a time axis. See
[ImageAxes](https://juliaimages.github.io/latest/imageaxes.html).

If you are passing the data via filenames, ensure that you have stored your images in an ITK-readable format. [NRRD](https://github.com/JuliaIO/NRRD.jl) is recommended. For those performing acquisition with Imagine, you can write out an [NRRD header](https://github.com/timholy/ImagineFormat.jl#converting-to-nrrd).

### Performing registration

There are two phases: writing out a deformation file and applying the deformation file (`warp`). You can combine the two with `register` by supplying an .


The process of registering images can be broken down into stages, and multiple stages can be cascaed
