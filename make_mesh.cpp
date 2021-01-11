/*
This program does the following in order:
1- Read a list of 3D points from input in a given text file.
   Each line in the input file is a point coordinates as: x y z
   (x,y) represents the pixel coordinates in the image I, where (0,0) is the 
   bottom-left point and (W-1,H-1) is the top-right point in I with
   height of H rows and width of W columns. The z coordinate is the gray-level
   value of the the pixel at location (x,y) in the image I. (i.e., I(x,y))
2- Generate the 2D Delaunay triangulation of the given points using their (x,y) coordinates and writes the 
   triangulation data into a given output file in 2D .OFF format as explained in https://en.wikipedia.org/wiki/OFF_(file_format)
   The best way to view the output OFF files is to use MeshLab software (https://www.meshlab.net/)  
3- Rasterize the mesh into an image as follows. For each face in the triangulation with vertices (x_i,y_i), (x_j,y_j), and (x_k,y_k), 
   form a planar interpolant that passes throught the points (x_i,y_i,z_i), (x_j,y_j,z_j), and (x_k,y_k,z_j).
   For each integer lattice point (x,y) in the mesh bounding box, find the face f containing (x,y) and
   use the planar interpolant for the face f to compute the function value I(x,y). Then, I(x,y) must be
   rounded to nearest integer to be used as the gray-level value in the reconstructed image. 
4- Write the sample values I(x,y) to the given output image file as the reconstructed image.
   If OpenCV is not installed, the output image is automatically stored in the .PNM format (https://en.wikipedia.org/wiki/Netpbm).
   The PNM image formats can be viewed online at: http://paulcuth.me.uk/netpbm-viewer/
   or, it can be converted to other image types at: https://convertio.co/
   If OpenCV is installed, then enable "#define OPENCV_INSTALLED" in file "include/config.hpp" to activate the
   OpenCV functionality. Then, you can use any common image formats (.png, .jpg, ...) in the output image file.
   
NOTES:
1- CGAL and GMP libraries are required for using this code.
   Since CGAL version 5.0, CGAL is header-only by default, which means that there is no need to compile CGAL or its libraries before it can be used.
   Only specify the include directory containing the "CGAL" folder in the compiler option if it's different from the system's include directory
   /usr/include in a Linux system.
2- OpenCV library is optional for generating the output images. Please refer to Step 4 above for further details.
3- CGAL::Projection_traits_xy_3 is used to be able to insert 3D (x,y,z) points into the 2D Delaunay triangulation.
   Normally, Delaunay triangulation only accepts 2D (x,y) points, but with the use of the CGAL::Projection_traits_xy_3
   class as the geometric traits for the trianulation, the 3D (x,y,z) points will be first projected to xy-plane and then 
   the projected 2D points (x,y) will be used for creating the triangulation. Also, the z values will be automatically 
   stored in the vertex data structure.

command line usage:
    ./make_mesh -i <input_points_file> -t <output_triangulation_file> -r <output_reconstructed_image>
For example:
    ./make_mesh -i input/lena_ed_4%.dat -t lena-tri.off -r lena-rec-img.pnm // if opencv not installed
    ./make_mesh -i input/lena_ed_4%.dat -t lena-tri.off -r lena-rec-img.png // if opencv installed
    
*/

#include <iostream>
#include <unistd.h>  //getopt
#include "include/config.hpp"
#include "include/MeshGeneration.hpp"

char* pointsFileName; //intput file for points list
char* triFileName;    //output file for storing triangulation mesh data in OFF format
char* imgFileName;    //output file for storing the reconstructed image

int main(int argc, char** argv)
{
    int opt;
    while ((opt = getopt(argc, argv, "i:t:r:")) != -1) 
    {
        switch (opt) 
        {
        case 'i':
            pointsFileName = optarg;
            break;
	    case 't':
	        triFileName = optarg;
	        break;
	    case 'r':
	        imgFileName = optarg;
	        break;
        }
    }

    //check is input/output files are specified
    if (pointsFileName == 0)
    {
        std::cerr << "ERROR: No input file for points list is specified!\n";
        std::cerr << "Usage: " << argv[0] << " -i <input_points_file> "
                                          << "-t <output_triangulation_file> "
                                          << "-r <output_reconstructed_image>\n";
        return -1;
    }
    else if (triFileName == 0 && imgFileName == 0 )
    {
        std::cerr << "ERROR: Please specify at least one output file name!\n";
        std::cerr << "Usage: " << argv[0] << " -i <input_points_file> "
                                          << "-t <output_triangulation_file> "
                                          << "-r <output_reconstructed_image>\n";
        return -1;
    }
   
    //default constructor to create the object
    MeshGeneration imgRec; 

    //read from input points list file and generate the triangulation mesh
    Triangulation tri = imgRec.genTri(pointsFileName); 

    //output triangulation file in OFF format
    if (triFileName != 0)
    {
        imgRec.writeTri(triFileName);
    }

    //output the reconstructed image
    if (imgFileName != 0)
    {
        #ifdef OPENCV_INSTALLED
        cv::Mat cvImg = imgRec.cvRasTri(tri);//rasterize triangulation to img data in 2D array
        imgRec.writeImg(imgFileName, cvImg); //write image in common formats as given by user

        #else
        std::vector<double> img = imgRec.rasTri(tri);//rasterize triangulation to img data in 1D array
        imgRec.writeImgPNM(imgFileName, img);       //write image in plain PGM format as defined by the Netpbm project
        #endif
    }

	return 0;
}
