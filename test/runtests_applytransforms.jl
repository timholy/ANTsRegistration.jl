using ANTsRegistration
using TestImages
using Images
using ImageView

#### Prepare test images
img = testimage("cameraman")
θ = deg2rad(15)
rot = imrotate(img, θ, center = true, axes(img))

# Pipeline
pipeline= [Stage(Global("Rigid"), MI(), (8,4,2), (3,2,1), (100,50,25), 1e-6, 10),
           Stage(Global("Affine"), MI(), (8,4,2), (3,2,1), (100,50,25), 1e-6, 10),
           Stage(SyN(), MI(), (4,2,1), (3,2,1), (30,30,20), 1e-6, 10)]

# Registration
tforms = register(img, rot, pipeline; seed = 1234)

# Transformation
tfms = [Tform(tforms[2]), Tform(tforms[1])]
imgw = applyTransforms(tfms, img, rot)
imshow(imgw)

# Inverse transformation
invtfms = [Tform(tforms[1], 1), Tform(tforms[3])]
imginv = applyTransforms(invtfms, img, imgw; verbose = true, suppressout = false) #FIXME
imshow(imginv)

# Apply transform to points
p = [Point(490,140,1,1), Point(259,407,1,1), Point(112, 173, 1, 1)]
pout = applyTransformsToPoints("test.csv", 2, invtfms, p)
