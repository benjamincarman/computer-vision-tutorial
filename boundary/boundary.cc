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
  Mat noise_free;
  Mat final;

  int low_h;
  int high_h;
  int low_s;
  int high_s;
  int low_v;
  int high_v;
};

bool isEdge(Mat &image, int i, int j)
{
  if (image.at<unsigned char>(i,j) == 255)
  {
    if (i > 0 && j > 0 && image.at<unsigned char>(i-1,j-1) == 0) return true;
    if (i > 0 && image.at<unsigned char>(i-1,j) == 0) return true;
    if (i > 0 && j < 255 && image.at<unsigned char>(i-1,j+1) == 0) return true;
    if (j > 0 && image.at<unsigned char>(i,j-1) == 0) return true;
    if (j < 255 && image.at<unsigned char>(i,j+1) == 0) return true;
    if (i < 255 && j > 0 && image.at<unsigned char>(i+1,j-1) == 0) return true;
    if (i < 255 && image.at<unsigned char>(i+1,j) == 0) return true;
    if (i < 255 && j < 255 && image.at<unsigned char>(i+1,j+1) == 0) return true;
  }

  if (image.at<unsigned char>(i,j) == 0)
  {
    if (i > 0 && j > 0 && image.at<unsigned char>(i-1,j-1) == 255) return true;
    if (i > 0 && image.at<unsigned char>(i-1,j) == 255) return true;
    if (i > 0 && j < 255 && image.at<unsigned char>(i-1,j+1) == 255) return true;
    if (j > 0 && image.at<unsigned char>(i,j-1) == 255) return true;
    if (j < 255 && image.at<unsigned char>(i,j+1) == 255) return true;
    if (i < 255 && j > 0 && image.at<unsigned char>(i+1,j-1) == 255) return true;
    if (i < 255 && image.at<unsigned char>(i+1,j) == 255) return true;
    if (i < 255 && j < 255 && image.at<unsigned char>(i+1,j+1) == 255) return true;
  }
  return false;
}

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
  namedWindow("Noise Free Image", WINDOW_AUTOSIZE);
  namedWindow("Final Image", WINDOW_AUTOSIZE);

  r.low_h = 0;
  r.high_h = 179;
  r.low_s = 159;//0;
  r.high_s = 255;//237;//170;
  r.low_v = 93;//1;
  r.high_v = 183;//255;

  inRange(r.hsv_img, Scalar(r.low_h, r.low_s, r.low_v), Scalar(r.high_h, r.high_s, r.high_v), r.hsv_img_thresholded);
  imshow("HSV Image Thresholded", r.hsv_img_thresholded);

  Mat labels, stats, centroids;
  int num_components = connectedComponentsWithStats(r.hsv_img_thresholded,
    labels, stats, centroids, 4);

  r.noise_free = Mat::zeros(r.hsv_img_thresholded.rows, r.hsv_img_thresholded.cols, CV_8U);
  for(int i = 0; i < r.noise_free.rows; i++) {
    for(int j = 0; j < r.noise_free.cols; j++) {
      int label = labels.at<int>(i,j);
      if (stats.at<int>(label, CC_STAT_AREA) > 2000 && stats.at<int>(label, CC_STAT_AREA) < 30000)
      {
        r.noise_free.at<unsigned char>(i,j) = 255;
      }
    }
  }

  r.final = Mat::zeros(r.noise_free.rows, r.noise_free.cols, CV_8U);

  //Get left boundary
  for(int i = 0; i < r.noise_free.rows; i++) {
    for(int j = 0; j < r.noise_free.cols; j++) {
      if (isEdge(r.noise_free, i, j))
      {
        r.final.at<unsigned char>(i, j) = 255;
        break;
      }
    }
  }

  //Get right boundary
  for(int i = 0; i < r.noise_free.rows; i++) {
    for(int j = r.noise_free.cols - 1; j >= 0; j--) {
      if (isEdge(r.noise_free, i, j))
      {
        r.final.at<unsigned char>(i, j) = 255;
        break;
      }
    }
  }

  //Get top boundary
  for(int j = 0; j < r.noise_free.cols; j++) {
    for(int i = 0; i < r.noise_free.rows; i++) {
      if (isEdge(r.noise_free, i, j))
      {
        r.final.at<unsigned char>(i, j) = 255;
        break;
      }
    }
  }

  //Get bottom boundary
  for(int j = 0; j < r.noise_free.cols; j++) {
    for(int i = r.noise_free.rows - 1; i >= 0; i--) {
      if (isEdge(r.noise_free, i, j))
      {
        r.final.at<unsigned char>(i, j) = 255;
        break;
      }
    }
  }

  imshow("Noise Free Image", r.noise_free);
  imshow("Final Image", r.final);

	waitKey(0);
  destroyWindow("Original Image");
  destroyWindow("HSV Image Thresholded");
  destroyWindow("Noise Free Image");
  destroyWindow("Final Image");

  return EXIT_SUCCESS;
}
