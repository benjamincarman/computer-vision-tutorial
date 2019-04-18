void get_lines(std::vector<Point> &boundary_points, Mat &image)
{
  Vec4f major;
  cv::fitLine(boundary_points, major, cv::DIST_L2, 1, 0.001, 0.001);

  std::cout << "Line[0]: " << major[0] << " Line[1]: " << major[1] << " Line[2]: " << major[2] << " Line[3]: " << major[3] << std::endl;
  /*
  Vec4f minor;
  minor[0] = -1 * major[1];
  minor[1] = major[0];
  minor[2] = major[2];
  minor[3] = major[3];
  */
  Point pmajor1, pmajor2, pminor1, pminor2;
  pmajor1.x = major[2] - 2000 * major[0];
  pmajor1.y = major[3] - 2000 * major[1];
  pmajor2.x = major[2] + 2000 * major[0];
  pmajor2.y = major[3] + 2000 * major[1];
  /*
  pminor1.x = minor[2] - 2000 * minor[0];
  pminor1.y = minor[3] - 2000 * minor[1];
  pminor2.x = minor[2] + 2000 * minor[0];
  pminor2.y = minor[2] + 2000 * minor[1];
  */
  cv::line(image, pmajor1, pmajor2, Scalar(255,0,0), 2);
  //cv::line(image, pminor1, pminor2, Scalar(100,0,0), 2);

  double slope = major[1]/major[0];
  double b = slope * -1 * major[2] + major[3];

  //vector<int> get_intersection_points(boundary_points, major[0], major[1]);

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



if (intersection1 < intersection2)
{
  std::vector<Point> set1 (boundary_points.begin() + intersection1, boundary_points.begin() + intersection2);

  std::vector<Point> set2 (boundary_points.begin(), boundary_points.begin() + intersection1);
  for (int i = intersection2; i < boundary_points.size(); i++)
  {
    set2.push_back(boundary_points[i]);
  }

  split_boundary.push_back(set1);
  split_boundary.push_back(set2);
}
else
{
  std::vector<Point> set1 (boundary_points.begin() + intersection2, boundary_points.begin() + intersection1);

  std::vector<Point> set2 (boundary_points.begin(), boundary_points.begin() + intersection2);
  for (int i = intersection1; i < boundary_points.size(); i++)
  {
    set2.push_back(boundary_points[i]);
  }

  split_boundary.push_back(set1);
  split_boundary.push_back(set2);
}






  /*
  std::list<std::pair<std::vector<Point>,Vec4f> > lines;

  for (int i = 0; i < boundary_lines.size(); i++)
  {
    Vec4f line;
    cv::fitLine(boundary_lines[i], line, cv::DIST_L2, 1, 0.001, 0.001);

    std::pair<std::vector<Point>,Vec4f> current (boundary_lines[i],line);
    lines.push_back(current);
  }
  std::cout << lines.size() << std::endl;

  int count = 0;
  for (std::list<std::pair<std::vector<Point>,Vec4f> >::iterator it = lines.begin(); it != lines.end(); it++)
  {
    for (std::list<std::pair<std::vector<Point>,Vec4f> >::iterator it2 = it++; it2 != lines.end(); it2++)
    {
      double vx1 = ((*it).second)[0];
      double vy1 = ((*it).second)[1];
      double ox1 = ((*it).second)[2];
      double oy1 = ((*it).second)[3];

      double vx2 = ((*it2).second)[0];
      double vy2 = ((*it2).second)[1];
      double ox2 = ((*it2).second)[2];
      double oy2 = ((*it2).second)[3];

      double dvx = abs(vx2 - vx1);
      double dvy = abs(vy2 - vy1);
      double dox = abs(ox2 - ox1);
      double doy = abs(oy2 - oy1);

      if (dvx < 0.25 && dvy < 0.25 && dox < 10 && doy < 10)
      {
        std::vector<Point> merged_vec = (*it).first;
        for (int i = 0; i < ((*it2).first).size(); i++)
        {
          merged_vec.push_back(((*it2).first)[i]);
        }

        Vec4f line;
        cv::fitLine(merged_vec, line, cv::DIST_L2, 1, 0.001, 0.001);

        std::pair<std::vector<Point>,Vec4f> merged_pair (merged_vec, line);

        it2 = lines.erase(it2);
        //it2--;
        //it = lines.erase(it);

        //it++;
        //lines.insert(it, merged_pair);
        //it--;
        //it = lines.erase(it);
        //it--;

        std::cout << "yo" << std::endl;

      }

    }
    */
