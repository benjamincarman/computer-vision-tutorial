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

void BGR_Slider(int, void* resource);
void HSV_Slider(int, void* resource);
void Mouse_Callback(int event, int x, int y, int flags, void* resource);

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

	namedWindow("Threshold Sliders", WINDOW_AUTOSIZE);
  namedWindow("Original Image", WINDOW_AUTOSIZE);
  imshow("Original Image", resource.bgr_img);
  namedWindow("BGR Thresholded Image", WINDOW_AUTOSIZE);
  namedWindow("HSV Thresholded Image", WINDOW_AUTOSIZE);

  resource.low_red = 255;
  resource.high_red = 0;
  resource.low_green = 255;
  resource.high_green = 0;
  resource.low_blue = 255;
  resource.high_blue = 0;

  resource.low_h = 179;
  resource.high_h = 0;
  resource.low_s = 255;
  resource.high_s = 0;
  resource.low_v = 255;
  resource.high_v = 0;

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

  //BGR_Slider(0, &resource);
  //HSV_Slider(0, &resource);

  setMouseCallback("Original Image", Mouse_Callback, &resource);

	waitKey(0);
	destroyWindow("Threshold Sliders");
  destroyWindow("Original Image");
  destroyWindow("BGR Thresholded Image");
  destroyWindow("HSV Thresholded Image");

  return EXIT_SUCCESS;
}

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

void Mouse_Callback(int event, int x, int y, int flags, void* resource)
{
  Threshold* r = (Threshold*)resource;

  if (event == EVENT_LBUTTONDOWN)
    {
      int min_b = 255;
      int max_b = 0;
      int min_g = 255;
      int max_g = 0;
      int min_r = 255;
      int max_r = 0;

      int min_h = 179;
      int max_h = 0;
      int min_s = 255;
      int max_s = 0;
      int min_v = 255;
      int max_v = 0;

      //Sample 10x10 array of pixels surrounding clicked location
      for (int i = x - 5; i < x + 5; i++)
      {
        for (int j = y - 5; j < y + 5; j++)
        {
          Vec3b bgr_point = r->bgr_img.at<Vec3b>(y, x);
          int blue = int(bgr_point.val[0]);
          int green = int(bgr_point.val[1]);
          int red = int(bgr_point.val[2]);
          std::cout << blue << " " << green << " " << red << std::endl;

          Vec3b hsv_point = r->hsv_img.at<Vec3b>(y, x);
          int hue = int(hsv_point.val[0]);
          int saturation = int(hsv_point.val[1]);
          int value = int(hsv_point.val[2]);
          std::cout << hue << " " << saturation << " " << value << std::endl;

          if (blue < min_b) {min_b = blue;}
          if (blue > max_b) {max_b = blue;}
          if (green < min_g) {min_g = green;}
          if (green > max_g) {max_g = green;}
          if (red < min_r) {min_r = red;}
          if (red > max_r) {max_r = red;}

          if (hue < min_h) {min_h = hue;}
          if (hue > max_h) {max_h = hue;}
          if (saturation < min_s) {min_s = saturation;}
          if (saturation > max_s) {max_s = saturation;}
          if (value < min_v) {min_v = value;}
          if (value > max_v) {max_v = value;}
        }
      }

      if (min_b < r->low_blue) {r->low_blue = min_b; setTrackbarPos("Low Blue Threshold Selector", "Threshold Sliders", r->low_blue);}
      if (max_b > r->high_blue) {r->high_blue = max_b; setTrackbarPos("High Blue Threshold Selector", "Threshold Sliders", r->high_blue);}
      if (min_g < r->low_green) {r->low_green = min_g; setTrackbarPos("Low Green Threshold Selector", "Threshold Sliders", r->low_green);}
      if (max_g > r->high_green) {r->high_green = max_g; setTrackbarPos("High Green Threshold Selector", "Threshold Sliders", r->high_green);}
      if (min_r < r->low_red) {r->low_red = min_r; setTrackbarPos("Low Red Threshold Selector", "Threshold Sliders", r->low_red);}
      if (max_r > r->high_red) {r->high_red = max_r; setTrackbarPos("High Red Threshold Selector", "Threshold Sliders", r->high_red);}

      if (min_h < r->low_h) {r->low_h = min_h; setTrackbarPos("Low Hue Threshold Selector", "Threshold Sliders", r->low_h);}
      if (max_h > r->high_h) {r->high_h = max_h; setTrackbarPos("High Hue Threshold Selector", "Threshold Sliders", r->high_h);}
      if (min_s < r->low_s) {r->low_s = min_s; setTrackbarPos("Low Saturation Threshold Selector", "Threshold Sliders", r->low_s);}
      if (max_s > r->high_s) {r->high_s = max_s; setTrackbarPos("High Saturation Threshold Selector", "Threshold Sliders", r->high_s);}
      if (min_v < r->low_v) {r->low_v = min_v; setTrackbarPos("Low Value Threshold Selector", "Threshold Sliders", r->low_v);}
      if (max_v > r->high_v) {r->high_v = max_v; setTrackbarPos("High Value Threshold Selector", "Threshold Sliders", r->high_v);}

      BGR_Slider(0, resource);
      HSV_Slider(0, resource);
    }
}
