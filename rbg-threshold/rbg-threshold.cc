#include "opencv2/imgproc.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui.hpp"
#include <cstdlib>
#include <iostream>

using namespace cv;

struct Threshold {
  Mat img;
  std::vector<Mat> planes;
  std::vector<Mat> planes_thresholded;
  Mat img_thresholded;
  int threshold_type;
  int low_red_threshold_value;
  int high_red_threshold_value;
  int low_green_threshold_value;
  int high_green_threshold_value;
  int low_blue_threshold_value;
  int high_blue_threshold_value;
};

void Slider(int, void* resource);

int main(int argc, char** argv)
{
  Threshold resource;

  //Read in image from input
  if (argc < 2)
  {
    std::cout << "Proper use is ./threshold <image-name> <opt-threshold-type>" << std::endl;
    return EXIT_FAILURE;
  }
	resource.img = imread(argv[1], -1);
	if (resource.img.empty()) {return -1;}

  split(resource.img, resource.planes);
  resource.planes_thresholded.push_back(resource.planes[0]);
  resource.planes_thresholded.push_back(resource.planes[1]);
  resource.planes_thresholded.push_back(resource.planes[2]);

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
  namedWindow("Image", WINDOW_AUTOSIZE);

  resource.low_red_threshold_value = 30;
  resource.high_red_threshold_value = 100;
  resource.low_green_threshold_value = 30;
  resource.high_green_threshold_value = 100;
  resource.low_blue_threshold_value = 30;
  resource.high_blue_threshold_value = 100;
  createTrackbar("Low Red Threshold Selector", "Threshold Sliders", &resource.low_red_threshold_value, 255, Slider, &resource);
  createTrackbar("High Red Threshold Selector", "Threshold Sliders", &resource.high_red_threshold_value, 255, Slider, &resource);
  createTrackbar("Low Green Threshold Selector", "Threshold Sliders", &resource.low_green_threshold_value, 255, Slider, &resource);
  createTrackbar("High Green Threshold Selector", "Threshold Sliders", &resource.high_green_threshold_value, 255, Slider, &resource);
  createTrackbar("Low Blue Threshold Selector", "Threshold Sliders", &resource.low_blue_threshold_value, 255, Slider, &resource);
  createTrackbar("High Blue Threshold Selector", "Threshold Sliders", &resource.high_blue_threshold_value, 255, Slider, &resource);

  Slider(0, &resource);

	waitKey(0);
	destroyWindow("Threshold Sliders");
  destroyWindow("Image");


  return EXIT_SUCCESS;
}

void Slider(int, void* resource)
{
  Threshold* r = (Threshold*)resource;
  inRange(r->img, Scalar(r->low_blue_threshold_value, r->low_green_threshold_value, r->low_red_threshold_value), Scalar(r->high_blue_threshold_value, r->high_green_threshold_value, r->high_red_threshold_value), r->img_thresholded);
  imshow("Image", r->img_thresholded);
}
