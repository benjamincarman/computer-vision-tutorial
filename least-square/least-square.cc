#include "opencv2/imgproc.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui.hpp"
#include <cstdlib>
#include <iostream>
#include <vector>
#include <list>

using namespace cv;

struct Config {
  Mat bgr_img;
  Mat hsv_img;
  Mat hsv_img_thresholded;
  Mat noise_free;
  Mat boundary;
  Mat thinned_boundary;
  Mat strip_tree;
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

double get_furthest_point(Vec4f major, const std::vector<Point> &line_points, int &index_of_furthest_point)
{
  double a = -1 * major[1]/major[0];
  double b = 1;
  double c = (major[1]/major[0]) * major[2] - major[3];

  double max_distance = 0;
  Point furthest_point;
  for (int i = 0; i < line_points.size(); i++)
  {
    double distance = abs(a*line_points[i].x + b*line_points[i].y + c)/sqrt(a*a + b*b);

    if (distance > max_distance)
    {
      max_distance = distance;
      furthest_point.x = line_points[i].x;
      furthest_point.y = line_points[i].y;
      index_of_furthest_point = i;
    }
  }

  return max_distance;
  //return furthest_point;
}

void split(std::vector<std::vector<Point > > &line_points)
{
  std::cout << "Split vec size: " << line_points.size() << std::endl;
  int last_index = line_points.size() - 1;
  Vec4f major;
  cv::fitLine(line_points[last_index], major, cv::DIST_L2, 1, 0.001, 0.001);
  int breaker;
  double distance = get_furthest_point(major, line_points[last_index], breaker);

  std::vector<Point> last = line_points[last_index];
  line_points.pop_back();
  std::vector<Point> first(last.begin(), last.begin() + breaker);
  std::vector<Point> second(last.begin() + breaker, last.end());
  line_points.push_back(first);
  if (distance < 5 && first.size() > 5)
  split(line_points);
  line_points.push_back(second);
  if (distance < 5 && second.size() > 5)
  split(line_points);
}

std::vector<std::vector<Point> > split_at_intersection_points(std::vector<Point> &boundary_points, Vec4f &major, Mat &image)
{
  std::vector<std::vector<Point> > split_boundary;

  //Get fitted line vector direction
  double vx = major[0];
  double vy = major[1];

  //Get origin fitted vector is respect to for line
  double ox = major[2];
  double oy = major[3];

  int intersection1;
  int intersection2;
  double min = 1000;
  double max = -1000;
  for (int i = 0; i < boundary_points.size(); i++)
  {
    //Get vector to current point from origin
    double px = boundary_points[i].x - ox;
    double py = boundary_points[i].y - oy;

    //Normalize vector
    double length = sqrt(px * px + py * py);
    px /= length;
    py /= length;

    //std::cout << "Point #: " << i << " x = " << thinned_boundary_points[i].x << "; y = " << thinned_boundary_points[i].y
    //          << "; Dot Product: " << x *  0.462352 + y * -0.886696 << std::endl;

    double dot_product = px *  vx + py * vy;
    if (dot_product < min)
    {
      intersection1 = i;
      min = dot_product;
    }
    if (dot_product > max)
    {
      intersection2 = i;
      max = dot_product;
    }
  }

  circle(image, boundary_points[intersection1], 5, Scalar(255, 0, 0));
  circle(image, boundary_points[intersection2], 5, Scalar(255, 0, 0));

  //Generate separate vectors of points and put into split_boundary
  if (intersection1 < intersection2)
  {
    std::vector<Point> set1 (boundary_points.begin() + intersection1, boundary_points.begin() + intersection2);

    std::vector<Point> set2 (boundary_points.begin() + intersection2, boundary_points.end());
    for (int i = 0; i < intersection1; i++)
    {
      set2.push_back(boundary_points[i]);
    }

    split_boundary.push_back(set1);
    split_boundary.push_back(set2);
  }
  else
  {
    std::vector<Point> set1 (boundary_points.begin() + intersection2, boundary_points.begin() + intersection1);

    std::vector<Point> set2 (boundary_points.begin() + intersection1, boundary_points.end());
    for (int i = 0; i < intersection2; i++)
    {
      set2.push_back(boundary_points[i]);
    }

    split_boundary.push_back(set1);
    split_boundary.push_back(set2);
  }

  //Debug info
  /*
  std::cout << "Size of split_boundary: " << split_boundary.size() << ", Size of boundary_points: " << boundary_points.size()
            << ", Size of split_boundary[0]: " << split_boundary[0].size() << ", Size of split_boundary[1]: " << split_boundary[1].size() << std::endl;
  std::cout << "Intersection1: " << intersection1 << ", Intersection2: " << intersection2 << std::endl;
  */

  return split_boundary;
}

void split_at_extreme_points(std::vector<std::vector<Point> > &boundary_lines, const std::vector<Point> &boundary_section, Mat &image)
{
  if (boundary_section.size() < 9)
  {
    boundary_lines.push_back(boundary_section);
    return;
  }

  Vec4f major;
/*
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  for (int i = 0; i < boundary_section.size(); i++)
  {
    std::cout << "Point " << i << ": x = " << boundary_section[i].x << ", y = " << boundary_section[i].y << std::endl;
    if (boundary_section.size() < 200) circle(image, boundary_section[i], 3, Scalar(255, 0, 0));
  }
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  */
  cv::fitLine(boundary_section, major, cv::DIST_L2, 1, 0.001, 0.001);

  Point pmajor1, pmajor2, pminor1, pminor2;
  pmajor1.x = major[2] - 2000 * major[0];
  pmajor1.y = major[3] - 2000 * major[1];
  pmajor2.x = major[2] + 2000 * major[0];
  pmajor2.y = major[3] + 2000 * major[1];

  //cv::line(image, pmajor1, pmajor2, Scalar(255,0,0), 2);

  //Get fitted line vector direction PERPENDICULAR
  double vx = major[1] * -1;
  double vy = major[0];

  //Get origin fitted vector is respect to for line
  double ox = major[2];
  double oy = major[3];

  std::list<std::pair<double,int> > dot_products;
  for (int i = 0; i < boundary_section.size(); i++)
  {
    //Get vector to current point from origin
    double px = boundary_section[i].x - ox;
    double py = boundary_section[i].y - oy;

    //Normalize vector
    double length = sqrt(px * px + py * py);
    px /= length;
    py /= length;

    //std::cout << "Point #: " << i << " x = " << thinned_boundary_points[i].x << "; y = " << thinned_boundary_points[i].y
    //          << "; Dot Product: " << x *  0.462352 + y * -0.886696 << std::endl;

    double dot_product = px *  vx + py * vy;
    std::pair<double,int> current_dot_prod (abs(dot_product),i);
    dot_products.push_back(current_dot_prod);
  }

  dot_products.sort();
  //circle(image, boundary_section[dot_products.back().second], 5, Scalar(255, 0, 0));
  //std::cout << "Biggest DP: " << dot_products.back().first << " Smallest DP: " << dot_products.front().first << std::endl;
  //circle(image, boundary_section[dot_products.front().second], 5, Scalar(255, 0, 0));

  //imshow("Strip Tree Lines", image);
  //std::cout << "Size: " << boundary_section.size() << std::endl;
  //waitKey();

  if (dot_products.back().first < 0.1)
  {
    boundary_lines.push_back(boundary_section);
    return;
  }
  else
  {
    int split_location = dot_products.back().second;
    std::vector<Point> set1 (boundary_section.begin(), boundary_section.begin() + split_location);
    std::vector<Point> set2 (boundary_section.begin() + split_location, boundary_section.end());

    split_at_extreme_points(boundary_lines, set1, image);
    split_at_extreme_points(boundary_lines, set2, image);
  }
}

void get_lines(std::vector<Point> &boundary_points, Mat &image)
{
  Vec4f major;
  cv::fitLine(boundary_points, major, cv::DIST_L2, 1, 0.001, 0.001);

  //Output debug info
  std::cout << "Line[0]: " << major[0] << " Line[1]: " << major[1] << " Line[2]: " << major[2] << " Line[3]: " << major[3] << std::endl;

  Point pmajor1, pmajor2, pminor1, pminor2;
  pmajor1.x = major[2] - 2000 * major[0];
  pmajor1.y = major[3] - 2000 * major[1];
  pmajor2.x = major[2] + 2000 * major[0];
  pmajor2.y = major[3] + 2000 * major[1];

  cv::line(image, pmajor1, pmajor2, Scalar(255,0,0), 2);

  std::vector<std::vector<Point> > split_boundary = split_at_intersection_points(boundary_points, major, image);

  //Do recursively
  //Pass a set of points, calculate the extreme
  //If extreme is below a certain threshold (or num points) push into a referenced vector
  //If not, create a set (vector) of points splitting at extremum, and continue recursively on each
  std::vector<std::vector<Point> > boundary_lines;
  split_at_extreme_points(boundary_lines, split_boundary[0], image);
  split_at_extreme_points(boundary_lines, split_boundary[1], image);

  for (int i = 0; i < boundary_lines.size(); i++)
  {
    Vec4f major;
    cv::fitLine(boundary_lines[i], major, cv::DIST_L2, 1, 0.001, 0.001);

    Point pmajor1, pmajor2, pminor1, pminor2;
    pmajor1.x = major[2] - 2000 * major[0];
    pmajor1.y = major[3] - 2000 * major[1];
    pmajor2.x = major[2] + 2000 * major[0];
    pmajor2.y = major[3] + 2000 * major[1];

    cv::line(image, pmajor1, pmajor2, Scalar(255,0,0), 1);
  }
  std::cout << "Size of boundary_lines: " << boundary_lines.size() << std::endl;
/*
  int index1 = -1;
  int index2 = -1;
  int count = 0;
  for (int i = 0; i < boundary_points.size(); i++)
  {
    double epsilon = 3.00;

    if (slope * boundary_points[i].x + b > boundary_points[i].y - epsilon &&
        slope * boundary_points[i].x + b < boundary_points[i].y + epsilon)
    {
        std::cout << "\nInitial Intersection:" << boundary_points[i] << std::endl;
        //circle(image, boundary_points[i], 5, Scalar(255, 0, 0));
        i += 15;
        count++;
        if (index1 == -1)
        {
          index1 = i;
        }
        else
        {
          index2 = i;
        }
    }
  }

  if (count != 2) {return;}

  std::vector<std::vector<Point> > line_points;

  std::vector<Point> first(boundary_points.begin(), boundary_points.begin() + index1);
  std::vector<Point> second(boundary_points.begin() + index1, boundary_points.begin() + index2);
  for (int i = index2; i < boundary_points.size(); i++)
  {
    first.push_back(boundary_points[i]);
  }
  line_points.push_back(first);
  //line_points.push_back(second);

  split(line_points);
  line_points.push_back(second);
  split(line_points);

  for (int i = 0; i < line_points.size(); i++)
  {
    Vec4f major;
    if (line_points[i].size() > 3)
    {
    cv::fitLine(line_points[i], major, cv::DIST_L2, 1, 0.001, 0.001);

    Point pmajor1, pmajor2;
    pmajor1.x = major[2] - 2000 * major[0];
    pmajor1.y = major[3] - 2000 * major[1];
    pmajor2.x = major[2] + 2000 * major[0];
    pmajor2.y = major[3] + 2000 * major[1];
    cv::line(image, pmajor1, pmajor2, Scalar(255,0,0), 2);
  }
  }

*/

  //std::cout << line_points.size() << std::endl;
  //cv::fitLine(line_points[0], major, cv::DIST_L2, 1, 0.001, 0.001);
  /*
  pmajor1.x = major[2] - 2000 * major[0];
  pmajor1.y = major[3] - 2000 * major[1];
  pmajor2.x = major[2] + 2000 * major[0];
  pmajor2.y = major[3] + 2000 * major[1];
  cv::line(image, pmajor1, pmajor2, Scalar(255,0,0), 2);
  */
  //Point breaker = get_furthest_point(major, line_points[0]);
  //circle(image, breaker, 5, Scalar(255, 0, 0));
  /*
  cv::fitLine(line_points[1], major, cv::DIST_L2, 1, 0.001, 0.001);
  pmajor1.x = major[2] - 2000 * major[0];
  pmajor1.y = major[3] - 2000 * major[1];
  pmajor2.x = major[2] + 2000 * major[0];
  pmajor2.y = major[3] + 2000 * major[1];
  cv::line(image, pmajor1, pmajor2, Scalar(255,0,0), 2);
  */
  //breaker = get_furthest_point(major, line_points[1]);
  //circle(image, breaker, 5, Scalar(255, 0, 0));
  //while (true)

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
    namedWindow("Strip Tree Lines", WINDOW_AUTOSIZE);
    namedWindow("Final Image", WINDOW_AUTOSIZE);
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
    /*
    int p1;
    int p2;
    double min = 1000;
    double max = -1000;
    for (int i = 0; i < thinned_boundary_points.size(); i++)
    {
      double x = thinned_boundary_points[i].x - 322.3;
      double y = thinned_boundary_points[i].y - 282.233;
      double length = sqrt(thinned_boundary_points[i].x * thinned_boundary_points[i].x + thinned_boundary_points[i].y * thinned_boundary_points[i].y);
      x /= length;
      y /= length;
      std::cout << "Point #: " << i << " x = " << thinned_boundary_points[i].x << "; y = " << thinned_boundary_points[i].y
                << "; Dot Product: " << x *  0.462352 + y * -0.886696 << std::endl;

      if (x *  0.462352 + y * -0.886696 < min)
      {
        p1 = i;
        min = x *  0.462352 + y * -0.886696;
      }
      if (x *  0.462352 + y * -0.886696 > max)
      {
        p2 = i;
        max = x *  0.462352 + y * -0.886696;
      }
    }
    */
    r.strip_tree = r.thinned_boundary.clone();
    get_lines(thinned_boundary_points, r.strip_tree);

    //circle(r.strip_tree, thinned_boundary_points[p1], 5, Scalar(255, 0, 0));
    //circle(r.strip_tree, thinned_boundary_points[p2], 5, Scalar(255, 0, 0));

    imshow("Noise Free Image", r.noise_free);
    imshow("Boundary", r.boundary);
    imshow("Thinned Boundary", r.thinned_boundary);
    imshow("Strip Tree Lines", r.strip_tree);



    r.line_length = 25;
    //createTrackbar("Minimum Line Length Selector", "Threshold Slider", &r.line_length, 100, slider, &r);
    //slider(0, &r);
    //Edge detection

    std::vector<Vec2f> lines;
    std::vector<Vec2f> unique_lines;
    HoughLines(r.boundary, lines, 1, 3*CV_PI/180, r.line_length);


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
        std::cout << rho << " " << theta << std::endl;
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

  	waitKey(0);
    destroyWindow("Original Image");
    destroyWindow("HSV Image Thresholded");
    destroyWindow("Noise Free Image");
    destroyWindow("Boundary");
    destroyWindow("Thinned Boundary");
    destroyWindow("Strip Tree Lines");
    destroyWindow("Final Image");
    //destroyWindow("Threshold Slider");

  }

  return EXIT_SUCCESS;
}
