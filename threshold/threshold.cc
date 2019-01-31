#include "opencv2/imgproc.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui.hpp"
#include <cstdlib>
#include <iostream>

using namespace cv;

struct Threshold {
  Mat img;
  Mat img_gray;
  Mat img_thresholded;
  int threshold_type;
  int threshold_value;
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

  //Make image grayscale
  cvtColor(resource.img, resource.img_gray, COLOR_BGR2GRAY);

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

	namedWindow("Threshold Slider", WINDOW_AUTOSIZE);
  namedWindow("Image", WINDOW_AUTOSIZE);

  resource.threshold_value = 0;
  createTrackbar("Threshold Selector", "Threshold Slider", &resource.threshold_value, 255, Slider, &resource);

  Slider(0, &resource);

	waitKey(0);
	destroyWindow("Threshold Slider");
  destroyWindow("Image");


  return EXIT_SUCCESS;
}

void Slider(int, void* resource)
{
  Threshold* r = (Threshold*)resource;
  threshold(r->img_gray, r->img_thresholded, r->threshold_value, 255, r->threshold_type);
  imshow("Image", r->img_thresholded);
}
