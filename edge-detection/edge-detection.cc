#include "opencv2/imgproc.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui.hpp"
#include <cstdlib>
#include <iostream>
#include <vector>

using namespace cv;

struct Threshold {
  Mat bgr_img;
  Mat hsv_img;
  Mat hsv_img_thresholded;
  Mat noise_free;
  Mat boundary;
  Mat corners;
  Mat final;

  int low_h;
  int high_h;
  int low_s;
  int high_s;
  int low_v;
  int high_v;

  int line_length;
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

void slider(int, void* resource)
{
  Threshold* r = (Threshold*)resource;

  //Edge detection
  std::vector<Vec2f> lines;
  HoughLines(r->boundary, lines, 1, CV_PI/180, r->line_length);

  r->final = r->boundary.clone();
  for (size_t i = 0; i < lines.size(); i++)
    {
      float rho = lines[i][0], theta = lines[i][1];
      Point pt1, pt2;
      double a = std::cos(theta), b = std::sin(theta);
      double x0 = a*rho, y0 = b*rho;
      pt1.x = cvRound(x0 + 1000*(-b));
      pt1.y = cvRound(y0 + 1000*(a));
      pt2.x = cvRound(x0 - 1000*(-b));
      pt2.y = cvRound(y0 - 1000*(a));
      line(r->final, pt1, pt2, Scalar(255,0,0), 2);
    }

  imshow("Final Image", r->final);

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
  namedWindow("Boundary", WINDOW_AUTOSIZE);
  namedWindow("Final Image", WINDOW_AUTOSIZE);
  namedWindow("Threshold Slider", WINDOW_AUTOSIZE);

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

  r.boundary = Mat::zeros(r.noise_free.rows, r.noise_free.cols, CV_8U);

  //Get left boundary
  for(int i = 0; i < r.noise_free.rows; i++) {
    for(int j = 0; j < r.noise_free.cols; j++) {
      if (isEdge(r.noise_free, i, j))
      {
        r.boundary.at<unsigned char>(i, j) = 255;
        break;
      }
    }
  }

  //Get right boundary
  for(int i = 0; i < r.noise_free.rows; i++) {
    for(int j = r.noise_free.cols - 1; j >= 0; j--) {
      if (isEdge(r.noise_free, i, j))
      {
        r.boundary.at<unsigned char>(i, j) = 255;
        break;
      }
    }
  }

  //Get top boundary
  for(int j = 0; j < r.noise_free.cols; j++) {
    for(int i = 0; i < r.noise_free.rows; i++) {
      if (isEdge(r.noise_free, i, j))
      {
        r.boundary.at<unsigned char>(i, j) = 255;
        break;
      }
    }
  }

  //Get bottom boundary
  for(int j = 0; j < r.noise_free.cols; j++) {
    for(int i = r.noise_free.rows - 1; i >= 0; i--) {
      if (isEdge(r.noise_free, i, j))
      {
        r.boundary.at<unsigned char>(i, j) = 255;
        break;
      }
    }
  }

  //Get edge lines
  r.line_length = 38;
  createTrackbar("Minimum Line Length Selector", "Threshold Slider", &r.line_length, 100, slider, &r);
  slider(0, &r);

  imshow("Noise Free Image", r.noise_free);
  imshow("Boundary", r.boundary);

	waitKey(0);
  destroyWindow("Original Image");
  destroyWindow("HSV Image Thresholded");
  destroyWindow("Noise Free Image");
  destroyWindow("Boundary");
  destroyWindow("Final Image");
  destroyWindow("Threshold Slider");

  return EXIT_SUCCESS;
}
