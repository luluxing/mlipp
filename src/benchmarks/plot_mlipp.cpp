#include <iostream>
#include <cfloat>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <chrono>

#include "mlipp_kd.h"
#include "point.h"
#include "multi_mlipp.h"

using namespace std;


int main(int argc, char** argv)
{
  if (argc != 2) {
    cout << "Usage: " << argv[0] << " <number_of_points>" << endl;
    return 1;
  }
  srand(9);

  int n = atoi(argv[1]);
  vector<Point<double>> data;
  for (int i = 0; i < n; ++i) {
    data.push_back((Point<double>){(double) rand() / RAND_MAX, (double) rand() / RAND_MAX});
  }

  MLIPP_KD<double> mlipp;
  mlipp.bulk_load(data.data(), data.size());
  mlipp.show();

  // MultiMlipp<double, true> multi_mlipp(2, Axis::X_AXIS, Partition::SPACE);
  // multi_mlipp.bulk_load(data.data(), data.size());
  // multi_mlipp.show();
  return 0;
}
