# Mesh Generation Project
## Summary ##
This program does the following in order:
1. Read a list of 3D points from the input text file. The points are assumed to be obtained as the result of a non-uniform sampling from a grayscale image. 
2. Generate the 2D Delaunay triangulation from the (x,y) coordinates of the input points.
3. Rasterize the mesh to integer grid points and write the reconstructed image to the output file.

For further implementation details please refer to the description included in file "make_mesh.cpp"

## Example ##
Consider the "boy.png" image from the "/images" folder as shown below: 

<img src="images/boy.png">

A list of non-uniformly-sampled 3D points from the "boy.png" image is provided in file "/input/boy_ed_4%.dat". The "4%" in the file name means that about 4% of total number of pixels (~= 0.04\*512\*512) in the image is used. Now, consider a command line as below:

./make_mesh -i /input/boy_ed_4%.dat -t output/boy-tri.off -r output/boy-rec-img.pnm

This command uses the "/input/boy_ed_4%.dat" file as the input to the "make_mesh" program and a Delaunay triangulation is generated as stored in "output/boy-tri.off" in .OFF format. A screenshot of the OFF file displyed by the MeshLab software is provided in file "output/boy-mesh.png" as shown below:

<img src="output/boy-mesh.png" width="512">

Then, the mesh is rasterized to integer grid points with the resolution same as the orignal image (i.e., 512\*512 grid points for "boy.png") to reconstruct the original image. The reconstructed image is written to file "output/boy-rec-img.pnm" in PNM format as shown below:

<img src="output/boy-rec-img.png">

The PNM image format can be viewed online at: http://paulcuth.me.uk/netpbm-viewer/ or can be converted to other image formats at: https://convertio.co/pnm-png/

If OpenCV is installed on your machine any common image formats can be used (i.e., JPG, PNG, ...) in the output image. To enable the OpenCV functionlity the "OPENVC_INSTALLED" macro must be first defined at the beginning of the file "make_mesh.cpp". More details are given inside the file "make_mesh.cpp".


