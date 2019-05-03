#include "opencv2/imgproc.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/calib3d/calib3d.hpp"
#include <cstdlib>
#include <iostream>
#include <vector>
#include <list>

#include "stop-detection.h"

using namespace cv;

struct Processed_Images {
  Mat bgr_img;
  Mat hsv_img;
  Mat hsv_img_thresholded;
  Mat noise_free;
  Mat boundary;
  Mat thinned_boundary;
  Mat boundary_corners;
  Mat perspective4_bgr;
  Mat perspective8_bgr;
  Mat perspective4_boundary;
  Mat perspective8_boundary;
  Mat perspective4_corners;
  Mat perspective8_corners;
  Mat strip_tree;
  Mat merged_strip_tree;
  Mat corners;
  Mat final;

  int selector;
};

//HSV Threshold Values
const int LOW_H = 0;
const int HIGH_H = 179;
const int LOW_S = 159;//0;
const int HIGH_S = 255;//237;//170;
const int LOW_V = 93;//1;
const int HIGH_V = 183;//255;

//Assume bgr_img is loaded
void generate_results(Processed_Images &results);

void display_results(Processed_Images &r);
void Image_Slider(int, void* results);
void Transformed_Image_Slider(int, void* results);

int main(int argc, char** argv)
{
  //Read in image from input
  if (argc < 2)
  {
    std::cout << "Proper use is ./boundary <image-name>" << std::endl;
    return EXIT_FAILURE;
  }

  for (int current_image = 1; current_image < argc; current_image++) {
    Processed_Images results;
    results.bgr_img = imread(argv[current_image], -1);
    if (results.bgr_img.empty()) {return -1;}

    generate_results(results);
    display_results(results);
  }

  return EXIT_SUCCESS;
}

void generate_results(Processed_Images &r)
{
  //Initialize HSV image
  cvtColor(r.bgr_img, r.hsv_img, COLOR_BGR2HSV);

  //Threshold for stop sign on HSV image
  inRange(r.hsv_img, Scalar(LOW_H, LOW_S, LOW_V), Scalar(HIGH_H, HIGH_S, HIGH_V), r.hsv_img_thresholded);

  //Get connected region
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


  //Get stop sign boundary
  r.boundary = Mat::zeros(r.noise_free.rows, r.noise_free.cols, CV_8U);

  std::vector<Point> boundary_points;

  //Find a Boundary
  bool done = false;
  for(int i = 0; i < r.noise_free.rows; i++) {
    if (done) break;
    for(int j = 0; j < r.noise_free.cols; j++) {
      if (isEdge(r.noise_free, i, j))
      {
        drawBoundary(r.noise_free, r.boundary, i, j, boundary_points);
        done = true;
        break;
      }
    }
  }

  //Generate thinned boundary
  r.thinned_boundary = r.boundary.clone();
  thin(r.thinned_boundary);

  Mat dummy = Mat::zeros(r.thinned_boundary.rows, r.thinned_boundary.cols, CV_8U);
  std::vector<Point> thinned_boundary_points;

  //Find thinned boundary points
  done = false;
  for(int i = 0; i < r.thinned_boundary.rows; i++) {
    if (done) break;
    for(int j = 0; j < r.thinned_boundary.cols; j++) {
      if (isEdge(r.thinned_boundary, i, j))
      {
        drawBoundary(r.thinned_boundary, dummy, i, j, thinned_boundary_points);
        done = true;
        break;
      }
    }
  }

  //Find boundary corners
  r.boundary_corners = r.thinned_boundary.clone();
  std::vector<Point> corner_points;
  get_corners(thinned_boundary_points, corner_points, r.boundary_corners, 0.09, 5);

  //Change perspective of bgr image using 4 points
  r.perspective4_bgr = r.bgr_img.clone();
  change_perspective4(corner_points, r.bgr_img, r.perspective4_bgr);

  //Change perspective of bgr image using 8 points
  r.perspective8_bgr = r.bgr_img.clone();
  change_perspective8(corner_points, r.bgr_img, r.perspective8_bgr);

  //Change perspective of boundary using 4 points
  r.perspective4_boundary = r.thinned_boundary.clone();
  change_perspective4(corner_points, r.thinned_boundary, r.perspective4_boundary);

  //Change perspective of boundary using all 8 points
  r.perspective8_boundary = r.thinned_boundary.clone();
  change_perspective8(corner_points, r.thinned_boundary, r.perspective8_boundary);

  r.perspective4_corners = r.perspective4_boundary.clone();
  r.perspective8_corners = r.perspective8_boundary.clone();

  get_perspective_corners(r.perspective4_corners);
  get_perspective_corners(r.perspective8_corners);

  //Generate boundary lines using strip trees
  r.strip_tree = r.thinned_boundary.clone();
  std::vector<std::vector<Point> > boundary_lines = get_lines(thinned_boundary_points, r.strip_tree);

  //Generate merged boundary lines
  r.merged_strip_tree = r.thinned_boundary.clone();
  merge_lines(boundary_lines, r.merged_strip_tree);

  //Generate boundary lines using Hough Lines
  int line_length = 25;

  std::vector<Vec2f> lines;
  std::vector<Vec2f> unique_lines;
  HoughLines(r.boundary, lines, 1, 3*CV_PI/180, line_length);

  int num_unique = 0;
  for (size_t i = 0; i < lines.size(); i++)
  {
    if (is_unique(lines[i], unique_lines))
    {
      unique_lines.push_back(lines[i]);
      num_unique++;
    }

    if (num_unique == 8) break;
  }

  r.final = r.boundary.clone();
  for (size_t i = 0; i < num_unique; i++)
    {
      float rho = unique_lines[i][0], theta = unique_lines[i][1];
      //std::cout << rho << " " << theta << std::endl;
      Point pt1, pt2;
      double a = std::cos(theta), b = std::sin(theta);
      double x0 = a*rho, y0 = b*rho;
      pt1.x = cvRound(x0 + 2000*(-b));
      pt1.y = cvRound(y0 + 2000*(a));
      pt2.x = cvRound(x0 - 2000*(-b));
      pt2.y = cvRound(y0 - 2000*(a));
      line(r.final, pt1, pt2, Scalar(255,0,0), 2);
    }
}

void display_results(Processed_Images &results)
{
  //namedWindow("Original Image", WINDOW_AUTOSIZE);
  //imshow("Original Image", results.bgr_img);

  namedWindow("Manipulated Images", WINDOW_AUTOSIZE);
  results.selector = 0;
  createTrackbar("Image Selector", "Manipulated Images", &results.selector, 8, Image_Slider, &results);
  imshow("Manipulated Images", results.bgr_img);

  namedWindow("Transformed Images", WINDOW_AUTOSIZE);
  results.selector = 0;
  createTrackbar("Image Selector", "Transformed Images", &results.selector, 5, Transformed_Image_Slider, &results);
  imshow("Transformed Images", results.perspective4_bgr);

  waitKey(0);
  //destroyWindow("Original Image");
  destroyWindow("Manipulated Images");
  destroyWindow("Transformed Images");
}

void Image_Slider(int, void* results)
{
  Processed_Images* r = (Processed_Images*)results;

  switch (r->selector)
  {
    case 0: imshow("Manipulated Images", r->bgr_img);
            break;
    case 1: imshow("Manipulated Images", r->hsv_img_thresholded);
            break;
    case 2: imshow("Manipulated Images", r->noise_free);
            break;
    case 3: imshow("Manipulated Images", r->boundary);
            break;
    case 4: imshow("Manipulated Images", r->thinned_boundary);
            break;
    case 5: imshow("Manipulated Images", r->boundary_corners);
            break;
    case 6: imshow("Manipulated Images", r->strip_tree);
            break;
    case 7: imshow("Manipulated Images", r->merged_strip_tree);
            break;
    case 8: imshow("Manipulated Images", r->final);
            break;
  }
}

void Transformed_Image_Slider(int, void* results)
{
  Processed_Images* r = (Processed_Images*)results;

  switch (r->selector)
  {
    case 0: imshow("Transformed Images", r->perspective4_bgr);
            break;
    case 1: imshow("Transformed Images", r->perspective8_bgr);
            break;
    case 2: imshow("Transformed Images", r->perspective4_boundary);
            break;
    case 3: imshow("Transformed Images", r->perspective8_boundary);
            break;
    case 4: imshow("Transformed Images", r->perspective4_corners);
            break;
    case 5: imshow("Transformed Images", r->perspective8_corners);
            break;
  }
}
