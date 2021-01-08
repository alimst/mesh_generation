//STL headers
#include <cmath>
#include <vector>
#include <string>
#include <cassert>
#include <fstream>
#include <map>

#ifdef OPENCV_INSTALLED
//OpenCV headers
#include <opencv2/core/utility.hpp>
#include <opencv2/core/mat.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#endif

//CGAL headers 
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Plane_3.h>

//Custom headers 
#include "rational.hpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Projection_traits_xy_3<Kernel> Gt;     //geometric traits for triangulation
typedef CGAL::Delaunay_triangulation_2<Gt> Triangulation;
typedef CGAL::Point_3<Kernel> Point;    // 3-D point
typedef CGAL::Bbox_2 Bbox;              // 2-D bounding box
typedef CGAL::Plane_3<Kernel> Plane;    // plane in 3-D space

typedef Triangulation::Finite_vertices_iterator fin_ver_iter;
typedef Triangulation::Finite_faces_iterator    fin_face_iter;
typedef Triangulation::Face Face;
typedef Triangulation::Face_handle Face_handle;
typedef Triangulation::Vertex_handle Vertex_handle;

//**********************************************************************
// Declarations
//**********************************************************************
class MeshGeneration
{
    //****************************
    //*** Private Data Members ***
    //****************************    
        Triangulation tri_;             // triangulation data
        std::vector<Point> pointsList_; // list of input points
        std::vector<double> img_;       // interpolant I(x,y) values of the reconstructed image
        Bbox B_;                    // triangulation bounding box
        int width_;                 // width of the image
        int height_;                // height of the image
        #ifdef OPENCV_INSTALLED
        cv::Mat cvImg_;             // 2D array for storing the reconstructed image
        #endif
        
    //********************************
    //*** Private Function Members ***
    //********************************
        // compute the planar interpolant of the given face.
        Plane computeFaceInterpolant ( Face_handle F);

        // return an integer flag representing the shape of that triangle face AND reorders the indicies of its 
        // vertices in a desired order such that the very bottom left vertex will always have the index 0 and the 
        // rest are 1 and 2 in counterclockwise order. So, face F is modified. 
        // There can be 4 face types:
        //      type 1: right-sided triangle
        //      type 2: left-sided triangle
        //      type 3: top two vertices have the same y
        //      type 4: bottom two vertices have the same y
        int getFaceType ( Face& F);

        // compute the bounding box of the given triangulation.
        Bbox computeBoundingBox() const;

        // compute the z value of the point (x,y,z) on the given plane.
        // lattice points (x,y) are integers.
        double computeZ ( const Plane& p, const int& x, const int& y) const;

    public:
        // read list of 3D points from standard input, build the triangulation, rasterize the mesh,
        // and write the reconstructed image to standard output
        MeshGeneration();

        // read points from input file and generate triangulation.
        Triangulation genTri(char* pointsFileName);

        // rasterize triangulation to find the I(x,y) values for the reconstructed image
        std::vector<double> rasTri(Triangulation& tri); //using 1D array without OpenCv

        // write the reconstructed image in PNM format (if opencv is not installed)
        void writeImgPNM( std::string outFile, std::vector<double>& img,
                           std::string type = "P2", int max_level = 255) const;

        //write triangulation data into file in OFF format
        void writeTri(char* triFileName) const;

        // return triangulation data
        Triangulation getTri() const;

        // return image data
        std::vector<double> getImg() const; //using 1D array

        // return image height
        int getHeight() const;

        // return image width
        int getWidth() const;

        #ifdef OPENCV_INSTALLED
        // rasterize triangulation to find the I(x,y) values for the reconstructed image
        cv::Mat cvRasTri(Triangulation& tri); //using 2D array from OpenCV

        //write the reconstructed image in any format (png, jpg, ...) using OpenCV library
        void writeImg(std::string outFile, cv::Mat& img) const;

        // return image data
        cv::Mat getCVImg() const; //using 2D array from opencv
        #endif

};

//**********************************************************************
// Definitions
//**********************************************************************
MeshGeneration::MeshGeneration()
{
    //data reset
    tri_.clear();
    pointsList_.clear();
    img_.clear();
    width_ = 0;
    height_= 0;
}

Triangulation MeshGeneration::getTri() const
{
    return tri_;
}

std::vector<double> MeshGeneration::getImg() const
{
    return img_;
}

#ifdef OPENCV_INSTALLED
cv::Mat MeshGeneration::getCVImg() const
{
    return cvImg_;
}

cv::Mat MeshGeneration::cvRasTri(Triangulation& tri)
{
    assert(tri.number_of_vertices() > 2);
    B_ = computeBoundingBox();// compute the bounding box B.
	width_ = B_.xmax() - B_.xmin() + 1.0;// width of the image (bounding box)
	height_ = B_.ymax() - B_.ymin() + 1.0;// height of the image (bounding box)

    //resize image data
    cvImg_ = cv::Mat::zeros(height_ , width_, CV_8UC1 );

	// Loop over only finite faces and compute interpolant over each face
	// NOTE: The interpolant values of the points inside the infinite faces were assigned to zero in 
    // initializing of img_.
    int faceCount = 0;
	for (fin_face_iter faceIter = tri.finite_faces_begin(); faceIter != tri.finite_faces_end(); ++faceIter)
	{
		Plane plane = computeFaceInterpolant( faceIter );// compute the planar interpolant
		int type = getFaceType ( *faceIter );// find the shape type of the face and reorder its indicies.
		
		// coordinates of the vertex with index 0:
		int x0 = faceIter->vertex(0)->point().x();
		int y0 = faceIter->vertex(0)->point().y();

		// coordinates of the vertex with index 1:
		int x1 = faceIter->vertex(1)->point().x();
		int y1 = faceIter->vertex(1)->point().y();

		// coordinates of the vertex with index 2:
		int x2 = faceIter->vertex(2)->point().x();
		int y2 = faceIter->vertex(2)->point().y();

		int height_1 = y1 - y0;
		int height_2 = y2 - y0;

		switch (type)
		{
        case 0: //invalid case - should never happen
            {
                assert(0);
                break;
            }
		case 1: // right-sided triangle
			{
            ++faceCount;
			rational slope_A = rational(x2-x0)/rational(y2-y0);// (Inverse of) the slope of the line between v0 and v2;
			rational slope_B = rational(x1-x0)/rational(y1-y0);// (Inverse of) the slope of the line between v0 and v1;
			rational slope_C = rational(x2-x1)/rational(y2-y1);// (Inverse of) the slope of the line between v1 and v2;

			// vertical scan from y0 to y2:
			for (int y_interval = 0; y_interval <= height_2; ++y_interval )
			{
				int y = y0 + y_interval;     // y-component of lattice point (x,y).				
				int start_x = rational::ceil(x0 + (y_interval*slope_A));
				int end_x;
				if ( y_interval <= height_1 )
				{
					end_x = rational::floor(x0 + (y_interval*slope_B));
				}else
				{
					end_x = rational::floor(x1 + ( (y_interval-height_1)*slope_C));
				}

				double z0 = computeZ(plane, start_x , y );
				int x_width = end_x - start_x;
				
				// horizontal scan from start_x to end_x:
				for ( int x_interval = 0; x_interval <= x_width; ++x_interval )
				{
					double z_interpolant = z0 - ( (x_interval)*(plane.a()/plane.c()) );
					int x = start_x + x_interval;// x-component of lattice point (x,y).
                    cvImg_.at<uchar>(B_.ymax() - y , x) = z_interpolant;
				}
	
			}
			break;
			}
			
		case 2: // left-sided triangle
			{ 
            ++faceCount;
			rational slope_A = rational(x2-x0)/rational(y2-y0);// (Inverse of) the slope of the line between v0 and v2;
			rational slope_B = rational(x1-x0)/rational(y1-y0);// (Inverse of) the slope of the line between v0 and v1;
			rational slope_C = rational(x2-x1)/rational(y2-y1);// (Inverse of) the slope of the line between v1 and v2;

			// vertical scan from y0 to y1:
			for (int y_interval = 0; y_interval <= height_1; ++y_interval )
			{
				int y = y0 + y_interval;     // y-component of lattice point (x,y).
				int end_x = rational::floor(x0 + (y_interval*slope_B));
				int start_x;
				if ( y_interval <= height_2 )
				{
					start_x = rational::ceil(x0 + (y_interval*slope_A));
				}else
				{
					start_x = rational::ceil(x2 + ( (y_interval-height_2)*slope_C));
				}

				double z0 = computeZ(plane, start_x , y );
				int x_width = end_x - start_x;
				
				// horizontal scan from start_x to end_x:
				for ( int x_interval = 0; x_interval <= x_width; ++x_interval )
				{
					double z_interpolant = z0 - ( x_interval*(plane.a()/plane.c()) );
					int x = start_x + x_interval;   // x-component of lattice point (x,y).
                    cvImg_.at<uchar>(B_.ymax() - y , x) = z_interpolant;
				}
			}
			break;
			}
		case 3: // (height_1 == height_2 , height_3 == 0)
			// top two vertices have the same y component.
			{
            ++faceCount;
			rational slope_A = rational(x2-x0)/rational(y2-y0);// (Inverse of) the slope of the line between v0 and v2;
			rational slope_B = rational(x1-x0)/rational(y1-y0);// (Inverse of) the slope of the line between v0 and v1;

			// vertical scan from y0 to y1==y2:
			for (int y_interval = 0; y_interval <= height_1; ++y_interval )
			{
				int y = y0 + y_interval;     // y-component of lattice point (x,y).
				int start_x = rational::ceil(x0 + (y_interval*slope_A));
				int end_x = rational::floor(x0 + (y_interval*slope_B));

				double z0 = computeZ(plane, start_x , y );
				int x_width = end_x - start_x;
				
				// horizontal scan from start_x to end_x:
				for ( int x_interval = 0; x_interval <= x_width; ++x_interval )
				{
					double z_interpolant = z0 - ( x_interval*(plane.a()/plane.c()) );
					int x = start_x + x_interval;// x-component of lattice point (x,y).
                    cvImg_.at<uchar>(B_.ymax() - y , x) = z_interpolant;

				}
			}
			break;
			}

		case 4: // ( height_1 == 0 , height_2 == height_3)
			// bottom two vertices have the same y component.
			{
            ++faceCount;
			rational slope_A = rational(x2-x0)/rational(y2-y0);// (Inverse of) the slope of the line between v0 and v2;
			rational slope_C = rational(x2-x1)/rational(y2-y1);// (Inverse of) the slope of the line between v1 and v2;

			// vertical scan from y0==y1 to y2:
			for (int y_interval = 0; y_interval <= height_2; ++y_interval )
			{
				int y = y0 + y_interval;     // y-component of lattice point (x,y).
				int start_x = rational::ceil(x0 + (y_interval*slope_A));
				int end_x = rational::floor(x1 + (y_interval*slope_C));

				double z0 = computeZ(plane, start_x , y );
				int x_width = end_x - start_x;
				
				// horizontal scan from start_x to end_x:
				for ( int x_interval = 0; x_interval <= x_width; ++x_interval )
				{
					double z_interpolant = z0 - ( x_interval*(plane.a()/plane.c()) );
					int x = start_x + x_interval;// x-component of lattice point (x,y).
                    cvImg_.at<uchar>(B_.ymax() - y , x) = z_interpolant;
                }
			}
			break;
			}
		}// end of switch	
	}// end of "face" for-loop 

    return cvImg_;
}

void MeshGeneration::writeImg(std::string outFile, cv::Mat& img) const
{
    cv::imwrite(outFile, img);
}
#endif

int MeshGeneration::getHeight() const
{
    return height_;
}

int MeshGeneration::getWidth() const
{
    return width_;
}

Triangulation MeshGeneration::genTri (char* pointsFileName)
{
    std::ifstream pointsFile(pointsFileName, std::ifstream::in);
	Point point;// 3-D point
	while (pointsFile >> point) // read points from file.
	{
        pointsList_.push_back(point);
		tri_.insert(point); // insert the point into triangulation
	}
    pointsFile.close();

    assert(pointsList_.size() != 0);
    assert(tri_.number_of_vertices() > 2);
    assert(pointsList_.size() == tri_.number_of_vertices() );

    return tri_;
}

std::vector<double> MeshGeneration::rasTri(Triangulation& tri)
{
    assert(tri.number_of_vertices() > 2);
    B_ = computeBoundingBox();// compute the bounding box B.
	width_ = B_.xmax() - B_.xmin() + 1.0;// width of the image (bounding box)
	height_ = B_.ymax() - B_.ymin() + 1.0;// height of the image (bounding box)

    //resize image data
    img_.resize(width_*height_ , -1.0);

	// Loop over only finite faces and compute interpolant over each face
	// NOTE: The interpolant values of the points inside the infinite faces were assigned to zero in 
    // initializing of img_.
    int faceCount = 0;
	for (fin_face_iter faceIter = tri.finite_faces_begin(); faceIter != tri.finite_faces_end(); ++faceIter)
	{
		Plane plane = computeFaceInterpolant( faceIter );// compute the planar interpolant
		int type = getFaceType ( *faceIter );// find the shape type of the face and reorder its indicies.
		
		// coordinates of the vertex with index 0:
		int x0 = faceIter->vertex(0)->point().x();
		int y0 = faceIter->vertex(0)->point().y();

		// coordinates of the vertex with index 1:
		int x1 = faceIter->vertex(1)->point().x();
		int y1 = faceIter->vertex(1)->point().y();

		// coordinates of the vertex with index 2:
		int x2 = faceIter->vertex(2)->point().x();
		int y2 = faceIter->vertex(2)->point().y();

		int height_1 = y1 - y0;
		int height_2 = y2 - y0;

		switch (type)
		{
        case 0: //invalid case - should never happen
            {
                assert(0);
                break;
            }
		case 1: // right-sided triangle
			{
            ++faceCount;
			rational slope_A = rational(x2-x0)/rational(y2-y0);// (Inverse of) the slope of the line between v0 and v2;
			rational slope_B = rational(x1-x0)/rational(y1-y0);// (Inverse of) the slope of the line between v0 and v1;
			rational slope_C = rational(x2-x1)/rational(y2-y1);// (Inverse of) the slope of the line between v1 and v2;

			// vertical scan from y0 to y2:
			for (int y_interval = 0; y_interval <= height_2; ++y_interval )
			{
				int y = y0 + y_interval;     // y-component of lattice point (x,y).				
				int start_x = rational::ceil(x0 + (y_interval*slope_A));
				int end_x;
				if ( y_interval <= height_1 )
				{
					end_x = rational::floor(x0 + (y_interval*slope_B));
				}else
				{
					end_x = rational::floor(x1 + ( (y_interval-height_1)*slope_C));
				}

				double z0 = computeZ(plane, start_x , y );
				int x_width = end_x - start_x;
				
				// horizontal scan from start_x to end_x:
				for ( int x_interval = 0; x_interval <= x_width; ++x_interval )
				{
					double z_interpolant = z0 - ( (x_interval)*(plane.a()/plane.c()) );
					int x = start_x + x_interval;// x-component of lattice point (x,y).
					img_[((y-B_.ymin())*width_)+(x-B_.xmin())] = z_interpolant; // corresponds to I(x,y).
				}
	
			}
			break;
			}
			
		case 2: // left-sided triangle
			{ 
            ++faceCount;
			rational slope_A = rational(x2-x0)/rational(y2-y0);// (Inverse of) the slope of the line between v0 and v2;
			rational slope_B = rational(x1-x0)/rational(y1-y0);// (Inverse of) the slope of the line between v0 and v1;
			rational slope_C = rational(x2-x1)/rational(y2-y1);// (Inverse of) the slope of the line between v1 and v2;

			// vertical scan from y0 to y1:
			for (int y_interval = 0; y_interval <= height_1; ++y_interval )
			{
				int y = y0 + y_interval;     // y-component of lattice point (x,y).
				int end_x = rational::floor(x0 + (y_interval*slope_B));
				int start_x;
				if ( y_interval <= height_2 )
				{
					start_x = rational::ceil(x0 + (y_interval*slope_A));
				}else
				{
					start_x = rational::ceil(x2 + ( (y_interval-height_2)*slope_C));
				}

				double z0 = computeZ(plane, start_x , y );
				int x_width = end_x - start_x;
				
				// horizontal scan from start_x to end_x:
				for ( int x_interval = 0; x_interval <= x_width; ++x_interval )
				{
					double z_interpolant = z0 - ( x_interval*(plane.a()/plane.c()) );
					int x = start_x + x_interval;   // x-component of lattice point (x,y).
					img_[((y-B_.ymin())*width_)+(x-B_.xmin())] = z_interpolant; // corresponds to I(x,y).
				}
			}
			break;
			}
		case 3: // (height_1 == height_2 , height_3 == 0)
			// top two vertices have the same y component.
			{
            ++faceCount;
			rational slope_A = rational(x2-x0)/rational(y2-y0);// (Inverse of) the slope of the line between v0 and v2;
			rational slope_B = rational(x1-x0)/rational(y1-y0);// (Inverse of) the slope of the line between v0 and v1;

			// vertical scan from y0 to y1==y2:
			for (int y_interval = 0; y_interval <= height_1; ++y_interval )
			{
				int y = y0 + y_interval;     // y-component of lattice point (x,y).
				int start_x = rational::ceil(x0 + (y_interval*slope_A));
				int end_x = rational::floor(x0 + (y_interval*slope_B));

				double z0 = computeZ(plane, start_x , y );
				int x_width = end_x - start_x;
				
				// horizontal scan from start_x to end_x:
				for ( int x_interval = 0; x_interval <= x_width; ++x_interval )
				{
					double z_interpolant = z0 - ( x_interval*(plane.a()/plane.c()) );
					int x = start_x + x_interval;// x-component of lattice point (x,y).
					img_[((y-B_.ymin())*width_)+(x-B_.xmin())] = z_interpolant; // corresponds to I(x,y).
				}
			}
			break;
			}

		case 4: // ( height_1 == 0 , height_2 == height_3)
			// bottom two vertices have the same y component.
			{
            ++faceCount;
			rational slope_A = rational(x2-x0)/rational(y2-y0);// (Inverse of) the slope of the line between v0 and v2;
			rational slope_C = rational(x2-x1)/rational(y2-y1);// (Inverse of) the slope of the line between v1 and v2;

			// vertical scan from y0==y1 to y2:
			for (int y_interval = 0; y_interval <= height_2; ++y_interval )
			{
				int y = y0 + y_interval;     // y-component of lattice point (x,y).
				int start_x = rational::ceil(x0 + (y_interval*slope_A));
				int end_x = rational::floor(x1 + (y_interval*slope_C));

				double z0 = computeZ(plane, start_x , y );
				int x_width = end_x - start_x;
				
				// horizontal scan from start_x to end_x:
				for ( int x_interval = 0; x_interval <= x_width; ++x_interval )
				{
					double z_interpolant = z0 - ( x_interval*(plane.a()/plane.c()) );
					int x = start_x + x_interval;// x-component of lattice point (x,y).
					img_[((y-B_.ymin())*width_)+(x-B_.xmin())] = z_interpolant; // corresponds to I(x,y).
                }
			}
			break;
			}
		}// end of switch	
	}// end of "face" for-loop 

    return img_;
}

Plane MeshGeneration::computeFaceInterpolant ( Face_handle F)
{
	Point p0 = F->vertex(0)->point();// point of the vertex with index 0 in a face.
	Point p1 = F->vertex(1)->point();// point of the vertex with index 1 in a face.
	Point p2 = F->vertex(2)->point();// point of the vertex with index 2 in a face.
	Plane plane(p0,p1,p2);           // construct a plane passing through p0, p1, and p2.
	return plane;
}

int MeshGeneration::getFaceType ( Face& F)
{
	int type = 0;	
	for (int i=0; i!=3; ++i)
	{
		int j = F.ccw(i); // returns i+1 in modulo 3
		int k = F.ccw(j); // returns j+1 in modulo 3

		// u, v, and w are three vertices of triangle in ccw order.
		Vertex_handle u = F.vertex( i );
 		Vertex_handle v = F.vertex( j );
		Vertex_handle w = F.vertex( k );

		// There are four cases based on the shape of the trinagle which are checked below:
		// check for case 1, case 2, and case 3 : 
		if ( ( u->point().y() < v->point().y() ) && ( u->point().y() < w->point().y() ) ) 
		{
			// case 1: right-sided triangle
			if ( ( v->point().y() < w->point().y() ))
			{
				F.set_vertex(0, u);// set the index of vertex u to 0.
				F.set_vertex(1, v);// set the index of vertex v to 1
				F.set_vertex(2, w);// set the index of vertex w to 2
				type = 1;
			}
			// case 2: left-sided triangle
			else if ( ( v->point().y() > w->point().y() ))
			{
				F.set_vertex(0, u);// set the index of vertex u to 0.
				F.set_vertex(1, v);// set the index of vertex v to 1
				F.set_vertex(2, w);// set the index of vertex w to 2
				type = 2;
			}
			// case 3: top two vertices have the same y component.
			else if ( ( v->point().y() == w->point().y() ))
			{
				F.set_vertex(0, u);// set the index of vertex u to 0.
				F.set_vertex(1, v);// set the index of vertex v to 1
				F.set_vertex(2, w);// set the index of vertex w to 2
				type = 3;
			}
		}
		// check for case 4: bottom two vertices have the same y component.
		else if ( ( u->point().y() == v->point().y() ))
		{
			F.set_vertex(0, u);// set the index of vertex u to 0.
			F.set_vertex(1, v);// set the index of vertex v to 1
			F.set_vertex(2, w);// set the index of vertex w to 2
			type = 4;
		}
	}
	return type;
}

Bbox MeshGeneration::computeBoundingBox() const
{
	// Loop over all finite vertices and compute the min and max coordinate components on x-y plane to form the bounding box.
	fin_ver_iter vertexIter = tri_.finite_vertices_begin();
	double minX = vertexIter->point().x();
	double minY = vertexIter->point().y();
	double maxX = vertexIter->point().x();
	double maxY = vertexIter->point().y(); 
	for (; vertexIter != tri_.finite_vertices_end(); ++vertexIter)
	{
		minX = std::min( minX , vertexIter->point().x() );
		maxX = std::max( maxX , vertexIter->point().x() );
		minY = std::min( minY , vertexIter->point().y() );
		maxY = std::max( maxY , vertexIter->point().y() );
	}
	Bbox B(minX,minY,maxX,maxY);
	return B;
}

double MeshGeneration::computeZ ( const Plane& P, const int& x, const int& y) const
{
	// General case of the plane equation is like ax+by+cz+d=0. 
    // Given x, y, and the plane, we can find the z component.
	double z = ((P.a())*x + (P.b())*y + P.d() )/(-(P.c()));// z=(ax+by+d)/(-c)
	return z;
} 

void MeshGeneration::writeImgPNM ( std::string outFile, std::vector<double>& img,
                                    std::string type, int max_level) const
{
    std::ofstream outImg (outFile , std::fstream::out);
    if ( !outImg.is_open() )
    {
        std::cerr << "Error! Failed to open file!\n";
        exit(1);
    }
    //PNM header
	outImg    << type << "\n"
		      << width_ << " "<< height_ << "\n"
		      << max_level << "\n";

    //pixels values
	for (int i=1; i <= height_; ++i) //height_
	{
		for (int j= 0; j!=width_; ++j) // width_
		{
			outImg << floor(img[(height_-i)*width_ + j]) <<" ";// pnm format only accepts integer numbers.
		}
		outImg<<"\n";
	}
    outImg.close();
}

void MeshGeneration::writeTri(char* triFileName) const
{
    std::ofstream triFile(triFileName , std::fstream::out);

    //OFF format header
    triFile << "OFF\n";
    triFile << tri_.number_of_vertices() << " " << tri_.number_of_faces() << " 0\n";
    triFile << '\n';

    //output vertices coordinates (with z=0) for 2D triangulation
    int verIndex=0;   //vertex index number
    std::map<Vertex_handle,int>  vertexMap; //vertex handle to index map
    vertexMap.clear();
    for(fin_ver_iter vertexIter=tri_.finite_vertices_begin(); vertexIter!=tri_.finite_vertices_end(); ++vertexIter)
    {
        vertexMap.insert(std::pair<Vertex_handle,int>(vertexIter,verIndex)); 
        triFile << vertexIter->point().x() << " " << vertexIter->point().y() << " 0\n";
        ++verIndex;
    }

    //output faces information
    triFile << '\n';
    for(fin_face_iter faceIter = tri_.finite_faces_begin(); faceIter != tri_.finite_faces_end(); ++faceIter)
    {
	    //find the three vertices of the face
        Vertex_handle v0 = faceIter->vertex(0);
        Vertex_handle v1 = faceIter->vertex(1);
        Vertex_handle v2 = faceIter->vertex(2);

        //find the indices of the vertices from the map container
        int ind0 = vertexMap.find(v0)->second;
        int ind1 = vertexMap.find(v1)->second;
        int ind2 = vertexMap.find(v2)->second;

        //output the indices of the three vertices v0,v1,v2 of the face
        triFile << "3 " << ind0 << " " << ind1 << " " << ind2 << '\n'; 
     }

     triFile.close();
}