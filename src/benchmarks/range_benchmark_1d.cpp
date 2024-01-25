#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <chrono>

#include "point.h"
#include "kdtree.h"
#include "lipp.h"

using namespace std;

uint compute_selectivity(const vector<Point>& data,
    const Point& query_min, const Point& query_max)
{
    uint cnt = 0;
    for (uint i = 0; i < data.size(); ++i)
    {
        if (data[i].x >= query_min.x 
            && data[i].x <= query_max.x)
            cnt++;
    }
    return cnt;
}

uint generate_query(const vector<Point>& data, float frac, 
    const Point& data_min, const Point& data_max,
    Point& query_min, Point& query_max)
{
    float dx = float(data_max.x - data_min.x);
    float qdx = frac * dx;
    int qx = rand() % int(dx - qdx);
    query_min = (Point){int(qx), int(qx)};
    query_max = (Point){int(qx + qdx), int(qx + qdx)};
    return compute_selectivity(data, query_min, query_max);
}

void run_lipp(vector<Point> data_pt, 
    vector<Point>& query_mins,
    vector<Point>& query_maxs,
    vector<uint>& counts)
{
    vector<pair<int, int>> data(data_pt.size());
    for (uint i = 0; i < data_pt.size(); ++i)
    {
        data[i] = make_pair(data_pt[i].x, data_pt[i].y);
    }

    LIPP<int, int> index;
    index.bulk_load(data.data(), data.size());

    for (uint i = 0; i < counts.size(); ++i)
    {
        vector<int> query_result(data.size());
        uint result_size = 0;
        auto start_time = std::chrono::high_resolution_clock::now();
        result_size = uint(index.range_query(query_result.data(), query_mins[i].x, query_maxs[i].x));
        auto end_time = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count() * 1e-9;
        if (result_size != counts[i])
            printf("Problem lipp, %d, %d, %d\n", i, int(result_size), counts[i]);
        else
            cout << "lipp," << counts[i] << "," << duration << endl;
    }
}

void run_kdtree(vector<Point> data, 
    vector<Point>& query_mins,
    vector<Point>& query_maxs,
    vector<uint>& counts)
{
    int max_depth = 0;
    KDNode* root = kd_bulk(data.data(), data.size(), 0, &max_depth);

    for (uint i = 0; i < counts.size(); ++i)
    {
        int count = 0;
        Point* query_result;
        auto start_time = std::chrono::high_resolution_clock::now();
        query_result = kd_range(root, query_mins[i], query_maxs[i], 
            0, max_depth, &count);
        auto end_time = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count() * 1e-9;
        free(query_result);
        if (uint(count) != counts[i])
            printf("Problem kdtree, %d, %d, %d\n", i, count, counts[i]);
        else
            cout << "kdtree," << counts[i] << "," << duration << endl;
    }
}

void run(uint num_points, uint num_fracs, uint iters_per_frac)
{
    Point data_min = (Point){RAND_MAX, RAND_MAX},
          data_max = (Point){0, 0};
    vector<Point> data_raw(num_points);
    for (uint i = 0; i < num_points; ++i)
    {
        data_raw[i].x = rand();
        data_raw[i].y = data_raw[i].x;
        if (data_raw[i].x < data_min.x)
            data_min.x = data_raw[i].x;
        if (data_raw[i].x > data_max.x)
            data_max.x = data_raw[i].x;
    }
    data_min.y = data_min.x;
    data_max.y = data_max.x;

    qsort(data_raw.data(), num_points, sizeof(Point), &compare_x);

    vector<int> idx;
    idx.push_back(0);
    for (uint i = 1; i < num_points; ++i)
        if (data_raw[i].x > data_raw[i - 1].x)
            idx.push_back(i);

    vector<Point> data;
    for (int i = 0; i < int(idx.size()); ++i)
        data.push_back(data_raw[idx[i]]);

    // random_shuffle(data.begin(), data.end());

    vector<Point> query_mins((num_fracs-1)*iters_per_frac), 
                  query_maxs((num_fracs-1)*iters_per_frac);
    vector<uint> counts((num_fracs-1)*iters_per_frac);
    int k = 0;
    for (uint i = 1; i < num_fracs; ++i)
    {
        for (uint j = 0; j < iters_per_frac; ++j)
        {
            counts[k] = generate_query(data, float(i) / float(num_fracs*100000),
                data_min, data_max, query_mins[k], query_maxs[k]);
            // printf("%d\n", counts[k]);
            k++;
        }
    }

    run_lipp(data, query_mins, query_maxs, counts);
    run_kdtree(data, query_mins, query_maxs, counts);

    return;
}



int main()
{
    run(1e7, 1e2, 1e2);
    return 0;
}