#include <iostream>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <time.h>
#include <chrono>

#include "kdtree.h"

using namespace std;

Point *
sequential_range_query(Point *points, int n, 
    Point min_point, Point max_point, int *count)
{
    int k = 0, size = 0;
    Point *result;
    bool *consistent = (bool *)malloc(sizeof(bool) * n);
    for (int i = 0; i < n; ++i)
    {
        if (points[i].x >= min_point.x
         && points[i].y >= min_point.y
         && points[i].x <= max_point.x
         && points[i].y <= max_point.y)
        {
            consistent[i] = true;
            size++;
        }
        else
            consistent[i] = false;

    }
    result = (Point *)malloc(sizeof(Point) * size);
    for (int i = 0; i < n; ++i)
    {
        if (consistent[i])
            result[k++] = points[i];
    }

    free(consistent);

    *count = size;
    return result;
}


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

    int max_depth = 0;
    root_bulk = kd_bulk(points, n, 0, &max_depth);

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

void
range_test(int n)
{
    KDNode *root;

    Point* points = (Point*)malloc(sizeof(Point) * n);
    for (int i = 0; i < n; ++i)
        points[i] = (Point){rand(), rand()};

    int max_depth = 0;
    root = kd_bulk(points, n, 0, &max_depth);

    Point min_point = (Point){rand(), rand()};
    Point max_point;
    do 
    {
        max_point = (Point){rand(), rand()};
    } 
    while (max_point.x <= min_point.x 
        || max_point.y <= min_point.y);

    int seq_count, kd_count;
    Point *seq_result = sequential_range_query(points, n, min_point, max_point, &seq_count);
    Point *kd_result;

    printf("Total size = %d, Result size = %d\n", n, seq_count);

    auto start_time = std::chrono::high_resolution_clock::now();

    kd_result = kd_range(root, min_point, max_point, max_depth, &kd_count);

    auto end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count() * 1e-9;
    cout << "Duration: range5 = " << duration << endl;

    if (seq_count == kd_count)
    {
        for (int i = 0; i < seq_count; ++i)
        {
            if (seq_result[i].x != kd_result[i].x 
                || seq_result[i].y != kd_result[i].y)
                printf("Not equal, %d\n", i);
        }
    }
    else
        printf("Range5: n_seq = %d, n_kd = %d\n", seq_count, kd_count);

    free(kd_result);
    free(seq_result);

}

// Define the main function
int main() {

    srand(time(NULL));

    /*run(1e6);
    for (int n = 5e6; n < 1e8; n += 5e6)
        run(n);*/

    range_test(1000000);

    return 0;
}
