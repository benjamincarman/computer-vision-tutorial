#include "opencv2/highgui/highgui.hpp"

using namespace cv;

int main(int argc, char** argv) 
{
	Mat img = imread(argv[1], -1);
	
	if (img.empty()) {return -1;}

	namedWindow("Example1", cv::WINDOW_AUTOSIZE);
	imshow("Example1", img);
	waitKey(0);
	destroyWindow("Example1");
}

/*
Compiled using the following:

gcc -v -std=c++11 -x c++ display.cc -I/usr/local/include/ -L/usr/lib/ -lstdc++ -L/usr/local/lib -lopencv_highgui -lopencv_core -lopencv_imgcodecs -o display


run with ./display logo.png

The parentheticals were added by me:

gcc -v (-std=c++11) (-x) (c++) display.cc -I/usr/local/include/ -L/usr/lib/ -lstdc++ -L/usr/local/lib -lopencv_highgui -lopencv_core (-lopencv_imgcodecs) -o display
*/
