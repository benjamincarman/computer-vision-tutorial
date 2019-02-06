#include "opencv2/imgproc.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui.hpp"
#include <cstdlib>
#include <iostream>

using namespace cv;

struct Threshold {
  Mat bgr_img;
  Mat bgr_img_thresholded;
  Mat hsv_img;
  Mat hsv_img_thresholded;

  int threshold_type;

  int low_red;
  int high_red;
  int low_green;
  int high_green;
  int low_blue;
  int high_blue;

  int low_h;
  int high_h;
  int low_s;
  int high_s;
  int low_v;
  int high_v;
};

void Type_Slider(int, void* resource);
void BGR_Slider(int, void* resource);
void HSV_Slider(int, void* resource);

int main(int argc, char** argv)
{
  Threshold resource;

  //Read in image from input
  if (argc < 2)
  {
    std::cout << "Proper use is ./threshold <image-name> <opt-threshold-type>" << std::endl;
    return EXIT_FAILURE;
  }
	resource.bgr_img = imread(argv[1], -1);
	if (resource.bgr_img.empty()) {return -1;}

  //Initialize HSV image
  cvtColor(resource.bgr_img, resource.hsv_img, COLOR_BGR2HSV);

  /* Default threshold type of 0 (binary) but allow other types be entered
     0: Binary
     1: Binary Inverted
     2: Threshold Truncated
     3: Threshold to Zero
     4: Threshold to Zero Inverted
  */
  resource.threshold_type = 0;
  if (argc > 2)
  {
    int input = atoi(argv[2]);
    if (input <= 4 & input >= 0) {resource.threshold_type = input;}
    else {
      std::cout << "Invalid threshold selctor entered" << std::endl;
      return EXIT_FAILURE;
    }
  }

	namedWindow("Threshold Sliders", WINDOW_AUTOSIZE);
  namedWindow("Original Image", WINDOW_AUTOSIZE);
  imshow("Original Image", resource.bgr_img);
  namedWindow("BGR Thresholded Image", WINDOW_AUTOSIZE);
  namedWindow("HSV Thresholded Image", WINDOW_AUTOSIZE);

  resource.low_red = 30;
  resource.high_red = 100;
  resource.low_green = 30;
  resource.high_green = 100;
  resource.low_blue = 30;
  resource.high_blue = 100;

  resource.low_h = 30;
  resource.high_h = 100;
  resource.low_s = 30;
  resource.high_s = 100;
  resource.low_v = 30;
  resource.high_v = 100;

  //createTrackbar("Threshold Type Selector", "Threshold Sliders", &resource.threshold_type, 4, Type_Slider, &resource);

  createTrackbar("Low Red Threshold Selector", "Threshold Sliders", &resource.low_red, 255, BGR_Slider, &resource);
  createTrackbar("High Red Threshold Selector", "Threshold Sliders", &resource.high_red, 255, BGR_Slider, &resource);
  createTrackbar("Low Green Threshold Selector", "Threshold Sliders", &resource.low_green, 255, BGR_Slider, &resource);
  createTrackbar("High Green Threshold Selector", "Threshold Sliders", &resource.high_green, 255, BGR_Slider, &resource);
  createTrackbar("Low Blue Threshold Selector", "Threshold Sliders", &resource.low_blue, 255, BGR_Slider, &resource);
  createTrackbar("High Blue Threshold Selector", "Threshold Sliders", &resource.high_blue, 255, BGR_Slider, &resource);

  createTrackbar("Low Hue Threshold Selector", "Threshold Sliders", &resource.low_h, 179, HSV_Slider, &resource);
  createTrackbar("High Hue Threshold Selector", "Threshold Sliders", &resource.high_h, 179, HSV_Slider, &resource);
  createTrackbar("Low Saturation Threshold Selector", "Threshold Sliders", &resource.low_s, 255, HSV_Slider, &resource);
  createTrackbar("High Saturation Threshold Selector", "Threshold Sliders", &resource.high_s, 255, HSV_Slider, &resource);
  createTrackbar("Low Value Threshold Selector", "Threshold Sliders", &resource.low_v, 255, HSV_Slider, &resource);
  createTrackbar("High Value Threshold Selector", "Threshold Sliders", &resource.high_v, 255, HSV_Slider, &resource);

  BGR_Slider(0, &resource);
  HSV_Slider(0, &resource);

	waitKey(0);
	destroyWindow("Threshold Sliders");
  destroyWindow("Original Image");
  destroyWindow("BGR Thresholded Image");
  destroyWindow("HSV Thresholded Image");

  return EXIT_SUCCESS;
}

/*
void Type_Slider(int, void* resource)
{
  Threshold* r = (Threshold*)resource;
  inRange(r->bgr_img, Scalar(r->low_blue, r->low_green, r->low_red), Scalar(r->high_blue, r->high_green, r->high_red), r->bgr_img_thresholded);
  imshow("BGR Thresholded Image", r->bgr_img_thresholded);
  inRange(r->hsv_img, Scalar(r->low_h, r->low_s, r->low_v), Scalar(r->high_h, r->high_s, r->high_v), r->hsv_img_thresholded);
  imshow("HSV Thresholded Image", r->hsv_img_thresholded);
}
*/

void BGR_Slider(int, void* resource)
{
  Threshold* r = (Threshold*)resource;

  if (r->high_blue < r->low_blue) {r->high_blue = r->low_blue; setTrackbarPos("High Blue Threshold Selector", "Threshold Sliders", r->high_blue);}
  if (r->high_green < r->low_green) {r->high_green = r->low_green; setTrackbarPos("High Green Threshold Selector", "Threshold Sliders", r->high_green);}
  if (r->high_red < r->low_red) {r->high_red = r->low_red; setTrackbarPos("High Red Threshold Selector", "Threshold Sliders", r->high_red);}

  inRange(r->bgr_img, Scalar(r->low_blue, r->low_green, r->low_red), Scalar(r->high_blue, r->high_green, r->high_red), r->bgr_img_thresholded);
  imshow("BGR Thresholded Image", r->bgr_img_thresholded);
}

void HSV_Slider(int, void* resource)
{
  Threshold* r = (Threshold*)resource;

  if (r->high_h < r->low_h) {r->high_h = r->low_h; setTrackbarPos("High Hue Threshold Selector", "Threshold Sliders", r->high_h);}
  if (r->high_s < r->low_s) {r->high_s = r->low_s; setTrackbarPos("High Saturation Threshold Selector", "Threshold Sliders", r->high_s);}
  if (r->high_v < r->low_v) {r->high_v = r->low_v; setTrackbarPos("High Value Threshold Selector", "Threshold Sliders", r->high_v);}

  inRange(r->hsv_img, Scalar(r->low_h, r->low_s, r->low_v), Scalar(r->high_h, r->high_s, r->high_v), r->hsv_img_thresholded);
  imshow("HSV Thresholded Image", r->hsv_img_thresholded);
}
