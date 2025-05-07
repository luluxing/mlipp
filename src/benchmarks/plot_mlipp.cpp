#include <iostream>
#include <cfloat>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <chrono>
#include <random>

#include "mlipp_kd.h"
#include "point.h"
#include "multi_mlipp.h"

using namespace std;

void gen_uniform_points(int n, std::vector<Point<double>>& points)
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dist_x(0.0, 1.0);
  std::uniform_real_distribution<> dist_y(0.0, 1.0);
  for (int i = 0; i < n; ++i) {
    double x = dist_x(gen);
    double y = dist_y(gen);
    // Ensure the points are within [0, 1] range
    x = std::max(0.0, std::min(1.0, x));
    y = std::max(0.0, std::min(1.0, y));
    points.push_back((Point<double>){x, y});
  }
}

void gen_normal_points(int n, std::vector<Point<double>>& points)
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<double> dist_x(0.5, 0.1); // Mean = 0.5, StdDev = 0.1
  std::normal_distribution<double> dist_y(0.5, 0.1); // Mean = 0.5, StdDev = 0.1

  for (int i = 0; i < n; ++i) {
    double x = dist_x(gen);
    double y = dist_y(gen);
    // Ensure the points are within [0, 1] range
    x = std::max(0.0, std::min(1.0, x));
    y = std::max(0.0, std::min(1.0, y));
    points.push_back((Point<double>){x, y});
  }
}

int main(int argc, char** argv)
{
  if (argc != 2) {
    cout << "Usage: " << argv[0] << " <number_of_points>" << endl;
    return 1;
  }
  srand(9);

  int n = atoi(argv[1]);
  vector<Point<double>> data;
  // gen_uniform_points(n, data);
  gen_normal_points(n, data);

  MLIPP_KD<double> mlipp;
  mlipp.bulk_load(data.data(), data.size());
  mlipp.show();

  // MultiMlipp<double, true> multi_mlipp(2, Axis::X_AXIS, Partition::SPACE);
  // multi_mlipp.bulk_load(data.data(), data.size());
  // multi_mlipp.show();
  return 0;
}
