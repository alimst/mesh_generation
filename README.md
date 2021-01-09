# Mesh Generation Project
## Summary ##
Thanks for your time to visit my Github page :)

Below is the summary of what the "make_mesh.cpp" program does in order:
1. Read a list of 3D points from the input text file. The points are assumed to be obtained as the result of a non-uniform sampling from a grayscale image. 
2. Generate the 2D Delaunay triangulation from the (x,y) coordinates of the input points.
3. Rasterize the mesh to integer grid points and write the reconstructed image to the output file.

For further implementation details please refer to the description included in file "make_mesh.cpp"

## Build Note ##
The "make_mesh.cpp" file can be built with g++ compiler either in Linux or with MinGW in Windows. The CGAL library is required for building the code (https://doc.cgal.org/latest/Manual/general_intro.html). The OpenCV library is optional for generating the output images of common formats. 

File ".vscodes/tasks.json" contains two build configurations: one for "windows" and one for "linux". The "windows" configuration shows an example on how to include required CGAL directories and link the libraries. The "linux" configuration provides an example on how to inlucde CGAL and OpenCV directories and link the required libraries. Locations of the include and library directories must be updated based on your own machine.

## Usage Example ##
Consider the "boy.png" image from the "/images" folder as shown below: 

<img src="images/boy.png">

A list of non-uniformly-sampled 3D points from the "boy.png" image is provided in file "/input/boy_ed_4%.dat". The "4%" in the file name means that about 4% of total number of pixels (~= 0.04\*512\*512) in the image is used. Now, consider a command line as below:

./make_mesh -i /input/boy_ed_4%.dat -t output/boy-tri.off -r output/boy-rec-img.pnm

This command uses the "/input/boy_ed_4%.dat" file as the input to the "make_mesh" program and a Delaunay triangulation is generated as stored in "output/boy-tri.off" in .OFF format. A screenshot of the OFF file displyed by the MeshLab software is provided in file "output/boy-mesh.png" as shown below:

<img src="output/boy-mesh.png" width="512">

Then, the mesh is rasterized to integer grid points with the resolution same as the orignal image (i.e., 512\*512 grid points for "boy.png") to reconstruct the original image. The reconstructed image is written to file "output/boy-rec-img.pnm" in PNM format as shown below:

<img src="output/boy-rec-img.png">

The PNM image format can be viewed online at: http://paulcuth.me.uk/netpbm-viewer/ or can be converted to other image formats at: https://convertio.co/pnm-png/. Just to grab your attention, this image was reconstructed using only 4% of all the pixels in the orignal image. 

If OpenCV is installed on your machine any common image formats can be used (i.e., JPG, PNG, ...) in the output image name. To enable the OpenCV functionality the "OPENVC_INSTALLED" macro must be first defined at the beginning of the file "make_mesh.cpp". More details are given inside the file "make_mesh.cpp".

## Contact ##
Please report any problems and bugs to ali.mostafavian@gmail.com
