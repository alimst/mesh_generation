# Makefile (keep this file in "make_mesh.cpp" directory)

# Set compiler command:
CXX = g++

# If OpenCV is installed set below to "YES", and update the include and lib directories paths down below, otherwise set to "NO".
OPENCV_INSTALLED = YES
ifeq ($(OPENCV_INSTALLED),YES) 
OPENCV_INC_PATH = /usr/local/include/opencv4
OPENCV_LIB_PATH = /usr/local/lib
OPENCV_LINKS= -lopencv_calib3d -lopencv_core -lopencv_features2d -lopencv_imgproc -lopencv_imgcodecs -lopencv_highgui
endif

# The C++ compiler and linker options:
# NOTE: Add include path for CGAL if different from default
# include directory /usr/include/
INCLUDES = -I $(OPENCV_INC_PATH)
LIBS = -L $(OPENCV_LIB_PATH)

CXXFLAGS = -O3 -std=c++17 $(INCLUDES) $(LIBS)
LINKS = -lgmp 

# Rules for building program:
make_mesh: make_mesh.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LINKS) $(OPENCV_LINKS)
	
# Rule for cleaning up files generated by compiling/linking:
clean:
	rm -f *.o make_mesh


