using ANTsRegistration
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
#tforms = register("tformout", img, rot, pipeline; seed = 1234) #save transform files in the hard drive

# Transformation
tfms = [Tform(tforms[2]), Tform(tforms[1])] #tforms[2]: warp, #tforms[1]: affine
imgw = applyTransforms(tfms, img, rot)

# Image Comparison
imgw1 = Gray{N0f8}.(imgw./typemax(UInt8))
imshow([img; imgw1])
imshow(RGB{N0f8}.(img, imgw1, zeros(Gray{N0f8}, size(img)))) #RGB

# Inverse transformation
invtfms = [Tform(tforms[1], 1), Tform(tforms[3])] #tforms[3]: invwarp, (tforms[1], 1): inv affine
imginv = applyTransforms(invtfms, img, imgw; verbose = true, suppressout = false) #FIXME
imshow(imginv)

# Apply transform to points
p = [Point(490,140,0,0), Point(259,407,0,0), Point(112, 173, 0, 0)] #Point(x, y, z, t) 
pout = applyTransformsToPoints(2, invtfms, p)
#pout = applyTransformsToPoints("fileout.csv", 2, invtfms, p) #save point coordinates in the hard drive
