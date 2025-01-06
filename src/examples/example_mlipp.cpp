#include <iostream>
#include <cfloat>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <chrono>

#include "mlipp_kd.h"
#include "point.h"

using namespace std;

template <typename T>
vector<Point<T>> seq_range_query(vector<Point<T>> points, 
    Point<T> min_key, Point<T> max_key)
{
    vector<Point<T>> result;
    for (uint i = 0; i < points.size(); ++i) {
        if (points[i].x >= min_key.x 
            && points[i].x <= max_key.x
            && points[i].y >= min_key.y
            && points[i].y <= max_key.y)
            result.push_back(points[i]);
    }
    return result;
}

template <typename T>
void run(int n)
{
    MLIPP_KD<T> mlipp;

    /* Generate random data */
    vector<Point<T>> data;
    
    double x1, x2, y1, y2;
    if (std::is_same<T, double>::value) {
      for (int i = 0; i < n; ++i) {
        data.push_back((Point<T>){(double) rand() / RAND_MAX, (double) rand() / RAND_MAX});
      }
      x1 = (double) rand() / RAND_MAX;
      x2 = (double) rand() / RAND_MAX;
      y1 = (double) rand() / RAND_MAX;
      y2 = (double) rand() / RAND_MAX;
  } else if (std::is_same<T, int>::value) {
      for (int i = 0; i < n; ++i) {
        data.push_back((Point<T>){rand(), rand()});
      }
      x1 = rand();
      x2 = rand();
      y1 = rand();
      y2 = rand();
    }

    /* Generate a random range */
    // double x1 = (double) rand() / RAND_MAX;
    // double x2 = (double) rand() / RAND_MAX;
    // double y1 = (double) rand() / RAND_MAX;
    // double y2 = (double) rand() / RAND_MAX;
    Point<T> min_point = (Point<T>){min(x1, x2), min(y1, y2)};
    Point<T> max_point = (Point<T>){max(x1, x2), max(y1, y2)};

    auto start_time = chrono::high_resolution_clock::now();

    /* Bulk loading the index */
    mlipp.bulk_load(data.data(), data.size());

    auto end_time = chrono::high_resolution_clock::now();
    auto construction_time = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count() * 1e-9;

    /* Compute index size */
    auto index_size = mlipp.index_size() / 1024.0 / 1024.0;

    start_time = chrono::high_resolution_clock::now();

    /* Point query */
    for (int i = 0; i < n; ++i)
        if (!mlipp.exists(data[i]))
            printf("Not found");

    end_time = chrono::high_resolution_clock::now();
    auto point_duration = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count() * 1e-9;

    start_time = chrono::high_resolution_clock::now();

    /* Range query */
    vector<Point<T>> seq_result = seq_range_query(data, min_point, max_point);
    vector<Point<T>> mlipp_result(seq_result.size());
    int result_size = 0;

    result_size = mlipp.range_query(min_point, max_point, mlipp_result.data());

    end_time = chrono::high_resolution_clock::now();
    auto range_duration = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count() * 1e-9;

    cout << "Construction time: " << construction_time << endl;
    cout << "Index size in MiB: " << index_size << endl;
    cout << "Point query duration: " << point_duration << endl;
    cout << "Range query duration: " << range_duration << endl;
    cout << "Range query duration / size: " << range_duration / result_size << endl;
}

int main()
{
    srand(9);
    run<int>(5e6);
    run<double>(5e6);
    return 0;
}
