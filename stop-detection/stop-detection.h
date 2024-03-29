using namespace cv;

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

  if (num_neighbors >= 1)
  {
    Mat labels;
    int num_components = connectedComponents(local_area, labels, 8);
    if (num_components == 2) //Only a foreground and background
    {
      return false;
    }
  }
  return true;

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
    //std::cout << count;

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

double get_difference(double angle1, double angle2)
{
  double difference1 = abs(angle1 - angle2);

  //If angles are in different quadrants
  if (angle1 * angle2 < 0)
  {
    double difference2;
    if (angle1 < 0)
    {
      difference2 = (CV_PI/2.0 + angle1) + (CV_PI/2.0 - angle2);
      return min(difference1, difference2);
    }
    else
    {
      difference2 = (CV_PI/2.0 + angle2) + (CV_PI/2.0 - angle1);
      return min(difference1, difference2);
    }
  }
  return difference1;
}

void get_corners(std::vector<Point> boundary_points, std::vector<Point> &corner_points, Mat &image, double threshold, int buffer)
{
  std::vector<double> slopes;

  //Get all out slopes from orgigin +- 6 pixels
  for (size_t i = 0; i < boundary_points.size(); i++)
  {
    if (i < buffer)
    {
      std::vector<Point> pts (boundary_points.begin(), boundary_points.begin() + (i + buffer));

      for (int j = boundary_points.size() - (buffer - i); j < boundary_points.size(); j++)
      {
        pts.push_back(boundary_points[j]);
      }

      Vec4f line; //use l2
      fitLine(pts, line, cv::DIST_L2, 1, 0.01, 0.01);

      double slope = atan2(line[1], line[0]);

      slopes.push_back(slope);
    }
    else if (i >= boundary_points.size() - buffer)
    {
      std::vector<Point> pts (boundary_points.begin() + (i - buffer), boundary_points.end());

      for (int j = 0; j < (buffer + 1) - (boundary_points.size() - i); j++)
      {
        pts.push_back(boundary_points[j]);
      }

      Vec4f line;
      fitLine(pts, line, cv::DIST_L2, 1, 0.01, 0.01);

      double slope = atan2(line[1], line[0]);

      slopes.push_back(slope);
    }
    else
    {
    std::vector<Point> pts (boundary_points.begin() + (i - buffer), boundary_points.begin() + (i + buffer));

    Vec4f line;
    fitLine(pts, line, cv::DIST_L2, 1, 0.01, 0.01);

    double slope = atan2(line[1], line[0]);

    slopes.push_back(slope);
    }
  }

  //Compute slope changes (gradient)
  std::vector<double> slope_change;

  for (int i = 0; i < slopes.size() - 1; i++)
  {
    double change = get_difference(slopes[i+1], slopes[i]);
    slope_change.push_back(change);
  }
  slope_change.push_back(get_difference(slopes[0], slopes[slopes.size() - 1]));

  std::vector<std::pair<Point,double> > corner_candidates;
  for (int i = 0; i < slope_change.size(); i++)
  {
    if (slope_change[i] >= threshold)
    {
      corner_candidates.push_back(std::make_pair(boundary_points[i], slope_change[i]));
    }
  }

  std::vector<bool> checked (corner_candidates.size(), false);

  for (int i = 0; i < corner_candidates.size(); i++)
  {
    if (checked[i]) continue;
    checked[i] = true;

    double max = corner_candidates[i].second;
    int index = i;

    for (int j = 0; j < corner_candidates.size(); j++)
    {
      if (checked[j]) continue;

      if (abs(corner_candidates[j].first.x - corner_candidates[i].first.x) < 15
          && abs(corner_candidates[j].first.y - corner_candidates[i].first.y) < 15)
      {
        checked[j] = true;

        if (corner_candidates[j].second > max)
        {
          max = corner_candidates[j].second;
          index = j;
        }
      }
    }

    corner_points.push_back(corner_candidates[index].first);
    circle(image, corner_candidates[index].first, 5, Scalar(255, 0, 0));
  }
}

void change_perspective4(std::vector<Point> corner_points, const Mat &source_image, Mat &dst_image)
{
  std::vector<Point2f> use_points;
  use_points.push_back(Point2f(corner_points[0]));
  use_points.push_back(Point2f(corner_points[1]));
  use_points.push_back(Point2f(corner_points[4]));
  use_points.push_back(Point2f(corner_points[5]));

  std::vector<Point2f> dst;

  dst.push_back(Point2f(139, 20));
  dst.push_back(Point2f(308, 20));
  dst.push_back(Point2f(308, 429));
  dst.push_back(Point2f(139, 429));

  Mat s(use_points);
  Mat d(dst);

  Mat transform = getPerspectiveTransform(s, d);
  warpPerspective(source_image, dst_image, transform, Size(447, 449));
}

void change_perspective8(std::vector<Point> corner_points, const Mat &source_image, Mat &dst_image)
{
  std::vector<Point2f> use_points;
  for (int i = 0; i < 8; i++)
  {
    use_points.push_back(Point2f(corner_points[i]));
  }

  std::vector<Point2f> dst;

  dst.push_back(Point2f(139, 20));
  dst.push_back(Point2f(308, 20));
  dst.push_back(Point2f(427, 140));
  dst.push_back(Point2f(427, 309));
  dst.push_back(Point2f(308, 429));
  dst.push_back(Point2f(139, 429));
  dst.push_back(Point2f(20, 309));
  dst.push_back(Point2f(20, 140));

  Mat s(use_points);
  Mat d(dst);

  Mat transform = findHomography(s, d);
  warpPerspective(source_image, dst_image, transform, Size(447, 449));
}

void get_perspective_corners(Mat &image)
{
  //Make boundary all same value
  for (int j = 0; j < image.cols; j++)
  {
    for (int i = 0; i < image.rows; i++)
    {
      if (image.at<unsigned char>(i,j) != 0)
      {
        image.at<unsigned char>(i,j) = 255;
      }
    }
  }

  //Generate thinned boundary
  thin(image);

  Mat dummy = Mat::zeros(image.rows, image.cols, CV_8U);
  std::vector<Point> bps;

  //Find thinned boundary points
  bool done = false;
  for(int i = 0; i < image.rows; i++) {
    if (done) break;
    for(int j = 0; j < image.cols; j++) {
      if (isEdge(image, i, j))
      {
        drawBoundary(image, dummy, i, j, bps);
        done = true;
        break;
      }
    }
  }

  //Find boundary corners
  std::vector<Point> corner_points;
  get_corners(bps, corner_points, image, 0.1, 7);
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

  //circle(image, boundary_points[intersection1], 5, Scalar(255, 0, 0));
  //circle(image, boundary_points[intersection2], 5, Scalar(255, 0, 0));

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

  if (dot_products.back().first < 0.8)
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

std::vector<std::vector<Point> > get_lines(std::vector<Point> &boundary_points, Mat &image)
{
  Vec4f major;
  cv::fitLine(boundary_points, major, cv::DIST_L2, 1, 0.001, 0.001);

  //Output debug info
  //std::cout << "Line[0]: " << major[0] << " Line[1]: " << major[1] << " Line[2]: " << major[2] << " Line[3]: " << major[3] << std::endl;
  /*
  Point pmajor1, pmajor2, pminor1, pminor2;
  pmajor1.x = major[2] - 2000 * major[0];
  pmajor1.y = major[3] - 2000 * major[1];
  pmajor2.x = major[2] + 2000 * major[0];
  pmajor2.y = major[3] + 2000 * major[1];

  cv::line(image, pmajor1, pmajor2, Scalar(255,0,0), 2);
  */
  std::vector<std::vector<Point> > split_boundary = split_at_intersection_points(boundary_points, major, image);

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

  return boundary_lines;
  //std::cout << "Size of boundary_lines: " << boundary_lines.size() << std::endl;
}

void merge_lines(std::vector<std::vector<Point> > &boundary_lines, Mat &image)
{
  for (int i = 0; i < boundary_lines.size() - 1; i++)
  {
    if (boundary_lines[i].size() == 0) continue;
    Vec4f line1;
    cv::fitLine(boundary_lines[i], line1, cv::DIST_L2, 1, 0.001, 0.001);
    double vx1 = line1[0];
    double vy1 = line1[1];
    double ox1 = line1[2];
    double oy1 = line1[3];

    for (int j = i + 1; j < boundary_lines.size(); j++)
    {
      if (boundary_lines[j].size() == 0) continue;
      Vec4f line2;
      cv::fitLine(boundary_lines[j], line2, cv::DIST_L2, 1, 0.001, 0.001);
      double vx2 = line2[0];
      double vy2 = line2[1];
      double ox2 = line2[2];
      double oy2 = line2[3];

      double dvx = abs(vx2 - vx1);
      double dvy = abs(vy2 - vy1);
      double dox = abs(ox2 - ox1);
      double doy = abs(oy2 - oy1);

      if (dvx < 0.5 && dvy < 0.5 && dox < 65 && doy < 70)
      {
        for (int k = 0; k < boundary_lines[j].size(); k++)
        {
          boundary_lines[i].push_back(boundary_lines[j][k]);
        }
        boundary_lines[j].clear();
      }
    }
  }

  for (int i = 0; i < boundary_lines.size(); i++)
  {
    if (boundary_lines[i].size() == 0) continue;

    Vec4f major;
    cv::fitLine(boundary_lines[i], major, cv::DIST_L2, 1, 0.001, 0.001);

    Point pmajor1, pmajor2, pminor1, pminor2;
    pmajor1.x = major[2] - 2000 * major[0];
    pmajor1.y = major[3] - 2000 * major[1];
    pmajor2.x = major[2] + 2000 * major[0];
    pmajor2.y = major[3] + 2000 * major[1];

    cv::line(image, pmajor1, pmajor2, Scalar(255,0,0), 1);
  }
}


//For HoughLines
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
