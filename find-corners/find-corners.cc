#include "opencv2/imgproc.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui.hpp"
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>

using namespace cv;

struct Config {
  Mat bgr_img;
  Mat hsv_img;
  Mat hsv_img_thresholded;
  Mat noise_free;
  Mat boundary;
  Mat thinned_boundary;
  Mat boundary_corners;
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
    if (i > 0 && j < image.cols && image.at<unsigned char>(i-1,j+1) == 0) return true;
    if (j > 0 && image.at<unsigned char>(i,j-1) == 0) return true;
    if (j < image.cols && image.at<unsigned char>(i,j+1) == 0) return true;
    if (i < image.rows && j > 0 && image.at<unsigned char>(i+1,j-1) == 0) return true;
    if (i < image.rows && image.at<unsigned char>(i+1,j) == 0) return true;
    if (i < image.rows && j < image.cols && image.at<unsigned char>(i+1,j+1) == 0) return true;
  }

  return false;
}

/*
could fit a quadraitc to area of 3 or 5 points. coefficiencts yield second partial derivatives for curvature

could use 4, 6, or 8 simple to thin the line and get rid of pixels then compute again

*/

void drawBoundary(Mat &source, Mat &target, int i, int j, std::vector<Point> &boundary_points)
{
  Point p;
  p.x = j;
  p.y = i;
  boundary_points.push_back(p);
  //std::cout<< i << " " << j << std::endl;
  target.at<unsigned char>(i,j) = 255;
  if (isEdge(source, i - 1, j) && target.at<unsigned char>(i - 1,j) == 0) drawBoundary(source, target, i - 1, j, boundary_points);
  if (isEdge(source, i - 1, j + 1) && target.at<unsigned char>(i - 1, j + 1) == 0) drawBoundary(source, target, i - 1, j + 1, boundary_points);
  if (isEdge(source, i, j + 1) && target.at<unsigned char>(i, j + 1) == 0) drawBoundary(source, target, i, j + 1, boundary_points);
  if (isEdge(source, i + 1, j + 1) && target.at<unsigned char>(i + 1, j + 1) == 0) drawBoundary(source, target, i + 1, j + 1, boundary_points);
  if (isEdge(source, i + 1, j) && target.at<unsigned char>(i + 1, j) == 0) drawBoundary(source, target, i + 1, j, boundary_points);
  if (isEdge(source, i + 1, j - 1) && target.at<unsigned char>(i + 1, j - 1) == 0) drawBoundary(source, target, i + 1, j - 1, boundary_points);
  if (isEdge(source, i, j - 1) && target.at<unsigned char>(i, j - 1) == 0) drawBoundary(source, target, i, j - 1, boundary_points);
  if (isEdge(source, i - 1, j - 1) && target.at<unsigned char>(i - 1, j - 1) == 0) drawBoundary(source, target, i - 1, j - 1, boundary_points);
}

bool changesConnectivity(Mat &image, int i, int j)
{
  int num_neighbors = 0;
  Mat local_area = Mat::zeros(3,3, CV_8U);

  for (int n = 0; n < 3; n++)
  {
    for (int m = 0; m < 3; m++)
    {
      //Keep the center a 0
      if (!(n == 1 && m == 1))
      {
        local_area.at<unsigned char>(n,m) = image.at<unsigned char>(i - 1 + n, j - 1 + m);

        if (local_area.at<unsigned char>(n,m) == 255)
        {
          num_neighbors++;
        }
      }
    }
  }
  /*
  if (local_area.at<unsigned char>(0,0) == 0 && local_area.at<unsigned char>(0,1) == 0 &&
      local_area.at<unsigned char>(0,2) == 0 && local_area.at<unsigned char>(1,0) == 0 &&
      local_area.at<unsigned char>(1,2) == 0 && local_area.at<unsigned char>(2,0) == 255 &&
      local_area.at<unsigned char>(2,1) == 255 && local_area.at<unsigned char>(2,2) == 255){
        std::cout << local_area << std::endl;
      }
  */
  if (num_neighbors > 1)
  {
    Mat labels;
    int num_components = connectedComponents(local_area, labels, 8);
    if (num_components == 2) //Only a foreground and background
    {
      return false;
    }
  }
  return true;
/*
  if (image.at<unsigned char>(i,j - 1) == 255 && image.at<unsigned char>(i,j + 1) == 255 &&
      image.at<unsigned char>(i - 1,j) == 0   && image.at<unsigned char>(i + 1,) == 0)
  {
    changes_connectivity = true;
  }

  if (image.at<unsigned char>(i,j - 1) == 255 && image.at<unsigned char>(i,j + 1) == 255 &&
      image.at<unsigned char>(i - 1,j) == 0   && image.at<unsigned char>(i + 1,) == 0)
  {
    changes_connectivity = true;
  }
*/
}

void thin(Mat &image)
{
  bool was_changed = true;
  int count = 0;
  while (was_changed)
  {
    was_changed = false;

    //Thin from top
    for (int j = 0; j < image.cols; j++)
    {
      for (int i = 0; i < image.rows; i++)
      {
        if (image.at<unsigned char>(i,j) == 255)
        {
          count++;
          if (!changesConnectivity(image, i, j))
          {
            image.at<unsigned char>(i,j) = 0;
            was_changed = true;
          }
          break;
        }
      }
    }
    std::cout << count;

    //Thin from right
    for (int i = 0; i < image.rows; i++)
    {
      for (int j = image.cols - 1; j >= 0; j--)
      {
        if (image.at<unsigned char>(i,j) == 255)
        {
          if (!changesConnectivity(image, i, j))
          {
            image.at<unsigned char>(i,j) = 0;
            was_changed = true;
          }
          break;
        }
      }
    }

    //Thin from bottom
    for (int j = 0; j < image.cols; j++)
    {
      for (int i = image.rows - 1; i >= 0; i--)
      {
        if (image.at<unsigned char>(i,j) == 255)
        {
          if (!changesConnectivity(image, i, j))
          {
            image.at<unsigned char>(i,j) = 0;
            was_changed = true;
          }
          break;
        }
      }
    }

    //Thin from left
    for (int i = 0; i < image.rows; i++)
    {
      for (int j = 0; j < image.cols; j++)
      {
        if (image.at<unsigned char>(i,j) == 255)
        {
          if (!changesConnectivity(image, i, j))
          {
            image.at<unsigned char>(i,j) = 0;
            was_changed = true;
          }
          break;
        }
      }
    }
  }
}

//atan2(delta y, delta x) 180 degree range in rad


void get_corners(std::vector<Point> boundary_points, std::vector<Point> &corner_points, Mat &image)
{

  std::vector<double> slopes;
  int buffer = 3;
  //Get all out slopes from orgigin +- 3 pixels
  for (size_t i = 0; i < boundary_points.size(); i++)
  {
    if (i < buffer)
    {
      std::vector<Point> pts (boundary_points.begin(), boundary_points.begin() + (i + buffer));

      for (int j = boundary_points.size() - (buffer - i); j < boundary_points.size(); j++)
      {
        pts.push_back(boundary_points[j]);
      }
      std::cout<< "#pts: " << pts.size()<< std::endl;
      Vec4f line; //use l2
      fitLine(pts, line, cv::DIST_L2, 1, 0.01, 0.01);
      /*
      Point pt1, pt2;
      pt1.x = line[2] - 100*line[0];
      pt1.y = line[3] - 100*line[1];
      pt2.x = line[2] + 100*line[0];
      pt2.y = line[3] + 100*line[1];
      cv::line(image, pt1, pt2, Scalar(255,0,0), 1);
      */
      double slope = atan(line[1]/line[0]) * 180 / CV_PI;// + CV_PI; //* 180 / CV_PI;//line[1]/line[0];
      //if(slope < 0) slope += 180;
      //slope *= (180 / CV_PI);
      slopes.push_back(slope);
      //std::cout << line[0] << " " << line[1] << " " << line[2] << " " << line[3] << " " << slope << std::endl;
    }
    else if (i >= boundary_points.size() - buffer)
    {
      std::vector<Point> pts (boundary_points.begin() + (i - buffer), boundary_points.end());

      for (int j = 0; j < (buffer + 1) - (boundary_points.size() - i); j++)
      {
        pts.push_back(boundary_points[j]);
      }
      std::cout<< "#pts: " << pts.size()<< std::endl;

      Vec4f line;
      fitLine(pts, line, cv::DIST_L2, 1, 0.01, 0.01);
      /*
      Point pt1, pt2;
      pt1.x = line[2] - 100*line[0];
      pt1.y = line[3] - 100*line[1];
      pt2.x = line[2] + 100*line[0];
      pt2.y = line[3] + 100*line[1];
      cv::line(image, pt1, pt2, Scalar(255,0,0), 1);
      */
      double slope = atan(line[1]/line[0]) * 180 / CV_PI;// + CV_PI; //* 180 / CV_PI;//line[1]/line[0];
      //if(slope < 0) slope += 180;
      //slope *= (180 / CV_PI);
      slopes.push_back(slope);
      //std::cout << line[0] << " " << line[1] << " " << line[2] << " " << line[3] << " " << slope << std::endl;
    }
    else
    {
    std::vector<Point> pts (boundary_points.begin() + (i - buffer), boundary_points.begin() + (i + buffer));
    std::cout<< "#pts: " << pts.size()<< std::endl;

    Vec4f line;
    fitLine(pts, line, cv::DIST_L2, 1, 0.01, 0.01);
    /*
    Point pt1, pt2;
    pt1.x = line[2] - 100*line[0];
    pt1.y = line[3] - 100*line[1];
    pt2.x = line[2] + 100*line[0];
    pt2.y = line[3] + 100*line[1];
    cv::line(image, pt1, pt2, Scalar(255,0,0), 1);
    */
    double slope = atan(line[1]/line[0]) * 180 / CV_PI;// + CV_PI; //* 180 / CV_PI;//line[1]/line[0];
    //if(slope < 0) slope += 180;
    //slope *= (180 / CV_PI);
    slopes.push_back(slope);
    //std::cout << line[0] << " " << line[1] << " " << line[2] << " " << line[3] << " " << slope << std::endl;
    }
    //std::cout << boundary_points[i] << std::endl;
  }

  //Compute slope changes (gradient)
  std::vector<double> slope_change;
  //Change over t=2
  /*
  slope_change.push_back((slopes[1] - slopes[slopes.size() - 1])/2);
  for (int i = 1; i < slopes.size() - 1; i++)
  {
    double change = (slopes[i+1] - slopes[i-1])/2;
    slope_change.push_back(change);
    std::cout << change << std::endl;
  }
  slope_change.push_back((slopes[0] - slopes[slopes.size() - 2])/2);
  */
  //Change over t=1
  for (int i = 0; i < slopes.size() - 1; i++)
  {
    double change = (slopes[i+1] - slopes[i]);
    slope_change.push_back(change);
    //std::cout << change << std::endl;
  }
  slope_change.push_back(slopes[0] - slopes[slopes.size() - 1]);

  int locs[8];
  for (int l = 0; l < 8; l++)
  {
    int max = 0;
    int loc;
    for (size_t i = 0; i < slope_change.size(); i++)
    {
      std::cout << "Iteration: " << l << ", i: " << i << ", Slope Change: " << slope_change[i] <<
                   ", Slope[i]: " << slopes[i] << ", Slope[i+1]: " << slopes[i+1] << std::endl;
      if (abs(slope_change[i]) > max)
      {
        max = abs(slope_change[i]);
        loc = i;
      }
    }

    std::cout << "Location " << loc << "was chosen" << std::endl;
    slope_change[loc] = 0;
    /*
    int start = (slope_change.size() + (loc - (slope_change.size()/10))) % slope_change.size();
    std::cout << "Loc: " << loc << " Start: " << start << std::endl;
    int count = 2 * (slope_change.size()/10);
    while (count > 0) {
      slope_change[start] = 0;

      count--;
      if (start + 1 >= slope_change.size())
      {
        start = 0;
      }
      else
      {
        start++;
      }
    }
    *////////
    /*
    for (int i = (slope_change.size() + (loc - (slope_change.size()/12))) % slope_change.size(); i < (loc + (slope_change.size()/12)) % slope_change.size(); i++)
    {
      slope_change[loc] = 0;
      std::cout << "Dude";
    }
    */
    int erasable_range = 1;//10;//slope_change.size()/9; //10
    if (loc < erasable_range)
    {
      for (int i = 0; i <= loc + erasable_range; i++)
      {
        slope_change[i] = 0;
      }

      for (int i = slope_change.size() - (erasable_range - loc); i < slope_change.size(); i++)
      {
        slope_change[i] = 0;
      }
    }
    else if (loc > slope_change.size() - (erasable_range + 1))
    {
      for (int i = loc - erasable_range; i <= slope_change.size(); i++)
      {
        slope_change[i] = 0;
      }

      for (int i = 0; i <= erasable_range - (slope_change.size() - loc); i++)
      {
        slope_change[i] = 0;
      }
    }
    else
    {
      for (int i = loc - erasable_range; i <= loc + erasable_range; i++)
      {
        slope_change[i] = 0;
      }
    }

    locs[l] = loc;
  }

  for (int i = 0; i < 8; i++)
  {
    std::cout << locs[i] << std::endl;
    corner_points.push_back(boundary_points[locs[i]]);

    circle(image, boundary_points[locs[i]], 5, Scalar(255, 0, 0));
  }



/*
  std::vector<Point> slope_time;
  for (int i = 0; i < slopes.size(); i++)
  {
    Point p;
    p.x = i * 1.0;
    p.y = cvRound(slopes[i]);
    slope_time.push_back(p);
    //std::cout << p.y << std::endl;
  }
  */
}


bool is_unique(Vec2f &line, std::vector<Vec2f> &unique_lines)
{
  float rho = line[0], theta = line[1];
  float rho_threshold = 66;
  float theta_threshold = 5 * (CV_PI / 180);

  for (size_t i = 0; i < unique_lines.size(); i++)
  {
    float rho_i = unique_lines[i][0], theta_i = unique_lines[i][1];
    if (abs(rho - rho_i) < rho_threshold && abs(theta - theta_i) < theta_threshold) return false;
  }

  return true;
}

int main(int argc, char** argv)
{
  //Read in image from input
  if (argc < 2)
  {
    std::cout << "Proper use is ./boundary <image-name>" << std::endl;
    return EXIT_FAILURE;
  }

  for (int i = 1; i < argc; i++) {
    Config r;
    r.bgr_img = imread(argv[i], -1);
  	if (r.bgr_img.empty()) {return -1;}

    //Initialize HSV image
    cvtColor(r.bgr_img, r.hsv_img, COLOR_BGR2HSV);

    namedWindow("Original Image", WINDOW_AUTOSIZE);
    imshow("Original Image", r.bgr_img);

    namedWindow("HSV Image Thresholded", WINDOW_AUTOSIZE);
    namedWindow("Noise Free Image", WINDOW_AUTOSIZE);
    namedWindow("Boundary", WINDOW_AUTOSIZE);
    namedWindow("Thinned Boundary", WINDOW_AUTOSIZE);
    namedWindow("Corner Points", WINDOW_AUTOSIZE);
    //namedWindow("Final Image", WINDOW_AUTOSIZE);
    //namedWindow("Threshold Slider", WINDOW_AUTOSIZE);

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

    imshow("Noise Free Image", r.noise_free);
    imshow("Boundary", r.boundary);

    r.thinned_boundary = r.boundary;
    thin(r.thinned_boundary);
    imshow("Thinned Boundary", r.thinned_boundary);


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


    r.boundary_corners = r.thinned_boundary;
    //Find corners
    std::vector<Point> corner_points;
    get_corners(thinned_boundary_points, corner_points, r.boundary_corners);

    imshow("Corner Points", r.boundary_corners);

/*
    r.line_length = 25;
    //createTrackbar("Minimum Line Length Selector", "Threshold Slider", &r.line_length, 100, slider, &r);
    //slider(0, &r);
    //Edge detection

    std::vector<Vec2f> lines;
    std::vector<Vec2f> unique_lines;
    HoughLines(r.thinned_boundary, lines, 1, 3*CV_PI/180, r.line_length);


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


    r.final = r.thinned_boundary.clone();
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

    imshow("Final Image", r.final);
*/
  	waitKey(0);
    destroyWindow("Original Image");
    destroyWindow("HSV Image Thresholded");
    destroyWindow("Noise Free Image");
    destroyWindow("Boundary");
    destroyWindow("Thinned Boundary");
    destroyWindow("Corner Points");

    //destroyWindow("Final Image");
    //destroyWindow("Threshold Slider");

  }

  return EXIT_SUCCESS;
}
