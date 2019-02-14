#include "opencv2/imgproc.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui.hpp"
#include <cstdlib>
#include <iostream>

using namespace cv;

struct Threshold {
  Mat bgr_img;
  Mat hsv_img;
  Mat hsv_img_thresholded;

  int low_h;
  int high_h;
  int low_s;
  int high_s;
  int low_v;
  int high_v;
};

int main(int argc, char** argv)
{
  Threshold r;

  //Read in image from input
  if (argc < 2)
  {
    std::cout << "Proper use is ./boundary <image-name>" << std::endl;
    return EXIT_FAILURE;
  }
	r.bgr_img = imread(argv[1], -1);
	if (r.bgr_img.empty()) {return -1;}

  //Initialize HSV image
  cvtColor(r.bgr_img, r.hsv_img, COLOR_BGR2HSV);

  namedWindow("Original Image", WINDOW_AUTOSIZE);
  imshow("Original Image", r.bgr_img);

  namedWindow("HSV Image Thresholded", WINDOW_AUTOSIZE);
  namedWindow("Final Image", WINDOW_AUTOSIZE);

  r.low_h = 0;
  r.high_h = 179;
  r.low_s = 0;
  r.high_s = 170;
  r.low_v = 1;
  r.high_v = 255;

  inRange(r.hsv_img, Scalar(r.low_h, r.low_s, r.low_v), Scalar(r.high_h, r.high_s, r.high_v), r.hsv_img_thresholded);
  imshow("HSV Image Thresholded", r.hsv_img_thresholded);

	waitKey(0);
  destroyWindow("Original Image");
  destroyWindow("HSV Image Thresholded");
  destroyWindow("Final Image");

  return EXIT_SUCCESS;
}
