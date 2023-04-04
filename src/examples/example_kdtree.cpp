#include <iostream>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <time.h>
#include <chrono>

#include "kdtree.h"

using namespace std;

void
run(int n)
{
    // Create a new kd-tree
    KDNode* root_insert = NULL;
    KDNode* root_bulk = NULL;

    Point* points = (Point*)malloc(sizeof(Point) * n);
    for (int i = 0; i < n; ++i)
        points[i] = (Point){rand(), rand()};

    auto start_time = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < n; ++i)
        root_insert = kd_insert(root_insert, points[i], 0);

    auto end_time = chrono::high_resolution_clock::now();
    auto duration_insert = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count() * 1e-9;

    start_time = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < n; ++i)
        kd_exists(root_insert, points[i], 0);

    end_time = chrono::high_resolution_clock::now();
    auto duration_scan_insert = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count() * 1e-9;

    start_time = std::chrono::high_resolution_clock::now();

    root_bulk = kd_bulk(points, n, 0);

    end_time = chrono::high_resolution_clock::now();
    auto duration_build = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count() * 1e-9;

    start_time = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < n; ++i)
        kd_exists(root_bulk, points[i], 0);

    end_time = chrono::high_resolution_clock::now();
    auto duration_scan_build = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count() * 1e-9;

    int max_depth_insert, max_depth_build;
    int sum_depth_insert, sum_depth_build;
    int sum_nodes_insert, sum_nodes_build;
    double avg_depth_insert, avg_depth_build;
    kd_depth(root_insert, 1, &max_depth_insert, &sum_depth_insert, &sum_nodes_insert);
    kd_depth(root_bulk, 1, &max_depth_build, &sum_depth_build, &sum_nodes_build);
    avg_depth_insert = double(sum_depth_insert) / double(sum_nodes_insert);
    avg_depth_build = double(sum_depth_build) / double(sum_nodes_build);
    cout << n << ", " 
         << duration_insert << ", " 
         << duration_scan_insert << ", " 
         << max_depth_insert << ", " 
         << avg_depth_insert << ", " 
         << duration_build << ", " 
         << duration_scan_build << ", " 
         << max_depth_build << ", " 
         << avg_depth_build 
    << endl;
}

// Define the main function
int main() {

    srand(time(NULL));

    run(1e6);
    for (int n = 5e6; n < 1e8; n += 5e6)
        run(n);

    return 0;
}
