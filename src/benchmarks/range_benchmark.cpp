#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <chrono>

#include "point.h"
#include "kdtree.h"
#include "mlipp_kd.h"
#include "mlipp_zorder.h"

using namespace std;

uint compute_selectivity(const vector<Point>& data,
    const Point& query_min, const Point& query_max)
{
    uint cnt = 0;
    for (uint i = 0; i < data.size(); ++i)
    {
        if (data[i].x >= query_min.x 
            && data[i].x <= query_max.x
            && data[i].y >= query_min.y
            && data[i].y <= query_max.y)
            cnt++;
    }
    return cnt;
}

uint generate_query(const vector<Point>& data, float frac, 
    const Point& data_min, const Point& data_max,
    Point& query_min, Point& query_max)
{
    float dx = float(data_max.x - data_min.x);
    float dy = float(data_max.y - data_min.y);
    float qdx = std::sqrt(frac) * dx;
    float qdy = std::sqrt(frac) * dx;
    int qx = rand() % int(dx - qdx);
    int qy = rand() % int(dy - qdy);
    query_min = (Point){int(qx), int(qy)};
    query_max = (Point){int(qx + qdx), int(qy + qdy)};
    return compute_selectivity(data, query_min, query_max);
}

void run_mlipp_kd(vector<Point> data, 
    vector<Point>& query_mins,
    vector<Point>& query_maxs,
    vector<uint>& counts)
{
    MLIPP_KD<int> index;
    index.bulk_load(data.data(), data.size());

    for (uint i = 0; i < counts.size(); ++i)
    {
        vector<Point> query_result;
        auto start_time = std::chrono::high_resolution_clock::now();
        index.rangeQuery(query_mins[i], query_maxs[i], query_result);
        auto end_time = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count() * 1e-9;
        if (query_result.size() != counts[i])
            printf("Problem mlipp_kd, %d, %d, %d\n", i, int(query_result.size()), counts[i]);
        else
            cout << "mlipp_kd," << counts[i] << "," << duration << endl;
    }
}

void run_mlipp_zorder(vector<Point> data, 
    vector<Point>& query_mins,
    vector<Point>& query_maxs,
    vector<uint>& counts)
{
    MLIPP_KD<int> index;
    index.bulk_load(data.data(), data.size());

    for (uint i = 0; i < counts.size(); ++i)
    {
        vector<Point> query_result(counts[i]);
        int result_size = 0;
        auto start_time = std::chrono::high_resolution_clock::now();
        result_size = index.range_query(query_mins[i], query_maxs[i], query_result.data());
        auto end_time = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count() * 1e-9;
        if (uint(result_size) != counts[i])
            printf("Problem mlipp_zorder, %d, %d, %d\n", i, result_size, counts[i]);
        else
            cout << "mlipp_zorder," << counts[i] << "," << duration << endl;
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
    Point data_min = (Point){(float) RAND_MAX, (float) RAND_MAX},
          data_max = (Point){0, 0};
    vector<Point> data(num_points);
    for (uint i = 0; i < num_points; ++i)
    {
        data[i].x = rand();
        data[i].y = rand();
        if (data[i].x < data_min.x)
            data_min.x = data[i].x;
        if (data[i].y < data_min.y)
            data_min.y = data[i].y;
        if (data[i].x > data_max.x)
            data_max.x = data[i].x;
        if (data[i].y > data_max.y)
            data_max.y = data[i].y;
    }

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

    run_mlipp_kd(data, query_mins, query_maxs, counts);
    run_mlipp_zorder(data, query_mins, query_maxs, counts);
    run_kdtree(data, query_mins, query_maxs, counts);

    return;
}



int main()
{
    run(1e7, 1e2, 1);
    return 0;
}