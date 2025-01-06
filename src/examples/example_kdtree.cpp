#include <iostream>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <time.h>
#include <chrono>

#include "kdtree.h"

using namespace std;

vector<Point> seq_range_query(vector<Point> points, 
    Point min_key, Point max_key)
{
    vector<Point> result;
    for (uint i = 0; i < points.size(); ++i) {
        if (points[i].x >= min_key.x 
            && points[i].x <= max_key.x
            && points[i].y >= min_key.y
            && points[i].y <= max_key.y)
            result.push_back(points[i]);
    }
    return result;
}

void
run(int n)
{
    KDNode* root = NULL;

    /* Generate random data */
    vector<Point> data;
    for (int i = 0; i < n; ++i)
      data.push_back((Point){(double) rand() / RAND_MAX, (double) rand() / RAND_MAX});

    /* Generate a random range */
    double x1 = (double) rand() / RAND_MAX;
    double x2 = (double) rand() / RAND_MAX;
    double y1 = (double) rand() / RAND_MAX;
    double y2 = (double) rand() / RAND_MAX;
    Point min_point = (Point){std::min(x1, x2), std::min(y1, y2)};
    Point max_point = (Point){std::max(x1, x2), std::max(y1, y2)};

    auto start_time = std::chrono::high_resolution_clock::now();

    /* Bulk loading the index */
    int max_depth = 0;
    root = kd_bulk(data.data(), data.size(), 0, &max_depth);

    auto end_time = chrono::high_resolution_clock::now();
    auto construction_time = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count() * 1e-9;

    /* Compute index size */
    auto index_size = kd_size(root) / 1024.0 / 1024.0;

    start_time = std::chrono::high_resolution_clock::now();

    /* Point query */
    for (int i = 0; i < n; ++i)
        if (!kd_exists(root, data[i], 0))
            printf("Not found");

    end_time = chrono::high_resolution_clock::now();
    auto point_duration = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count() * 1e-9;

    start_time = std::chrono::high_resolution_clock::now();

    int kd_count;
    vector<Point> seq_result = seq_range_query(data, min_point, max_point);
    // Point *kd_result;

    // kd_result = kd_range(root, min_point, max_point, 0, max_depth, &kd_count);

    end_time = chrono::high_resolution_clock::now();
    auto range_duration = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count() * 1e-9;

    cout << "Construction time: " << construction_time << endl;
    cout << "Index size in MiB: " << index_size << endl;
    cout << "Point query duration: " << point_duration << endl;
    cout << "Range query duration: " << range_duration << endl;
    cout << "Range query duration / size: " << range_duration / kd_count << endl;
}

int main()
{
    srand(9);
    run(5e6);
    return 0;
}
