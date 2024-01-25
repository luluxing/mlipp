#include <iostream>
#include <cfloat>
#include <stdlib.h>
#include <time.h>
#include <chrono>

#include "mlipp_kd.h"
#include "point.h"

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

int seq_knn_query(vector<Point> points, Point query_point, int k, 
    Point* results, double* distances)
{
    double maxDist = DBL_MAX;
    int n = 0;
    for (uint j = 0; j < points.size(); ++j) {
        Point p = points[j];
        double dist = pow(p.x - query_point.x, 2) 
                    + pow(p.y - query_point.y, 2);
        if (dist < maxDist) {
            int i = n;
            while (i > 0 && dist < distances[i-1])
                i--;
            int shift_count = n < k ? n - i : n - i - 1;
            memmove(&results[i + 1], &results[i], sizeof(Point)*shift_count);
            memmove(&distances[i + 1], &distances[i], sizeof(double)*shift_count);
            results[i] = p;
            distances[i] = dist;
            if (n < k)
                n++;
            if (n == k)
                maxDist = distances[n - 1];
        }
    }
    for (int i = 0; i < n; ++i)
        distances[i] = sqrt(distances[i]);
    return n;
}

void
run(int n)
{
    MLIPP_KD<int> lipp_insert;
    MLIPP_KD<int> lipp_bulk;

    vector<Point> data;
    for (int i = 0; i < n; ++i)
        data.push_back((Point){rand(), rand()});

    auto start_time = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < n; ++i)
        lipp_insert.insert(data[i]);

    auto end_time = chrono::high_resolution_clock::now();
    auto duration_insert = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count() * 1e-9;

    start_time = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < n; ++i)
        if (!lipp_insert.exists(data[i]))
            printf("Not found");

    end_time = chrono::high_resolution_clock::now();
    auto duration_scan_insert = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count() * 1e-9;

    start_time = std::chrono::high_resolution_clock::now();

    lipp_bulk.bulk_load(data.data(), data.size());

    end_time = chrono::high_resolution_clock::now();
    auto duration_build = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count() * 1e-9;

    start_time = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < n; ++i)
        if (!lipp_bulk.exists(data[i]))
            printf("Not found");

    end_time = chrono::high_resolution_clock::now();
    auto duration_scan_build = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count() * 1e-9;

    int max_depth_insert, max_depth_build;
    double avg_depth_insert, avg_depth_build;
    lipp_insert.print_depth(&max_depth_insert, &avg_depth_insert);
    lipp_bulk.print_depth(&max_depth_build, &avg_depth_build);
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
    MLIPP_KD<int> mlipp;

    vector<Point> data;
    for (int i = 0; i < n; ++i)
        data.push_back((Point){rand(), rand()});

    mlipp.bulk_load(data.data(), data.size());

    int x1 = rand();
    int x2 = rand();
    int y1 = rand();
    int y2 = rand();
    Point min_point = (Point){std::min(x1, x2), std::min(y1, y2)};
    Point max_point = (Point){std::max(x1, x2), std::max(y1, y2)};

    std::vector<Point> seq_results = seq_range_query(data, min_point, max_point);
    std::vector<Point> mlipp_results(seq_results.size());
    int result_size = 0;

    printf("Total size = %d, Result size = %d\n", n, int(seq_results.size()));

    auto start_time = std::chrono::high_resolution_clock::now();

    result_size = mlipp.range_query(min_point, max_point, mlipp_results.data());

    auto end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count() * 1e-9;
    cout << "Duration: range = " << duration << endl;
    cout << "Duration / size: range = " << duration / result_size << endl;


    qsort(seq_results.data(), seq_results.size(), sizeof(Point), &compare_pt);
    qsort(mlipp_results.data(), result_size, sizeof(Point), &compare_pt);

    if (int(seq_results.size()) == result_size)
    {
        for (uint i = 0; i < seq_results.size(); ++i)
        {
            if (seq_results[i].x != mlipp_results[i].x 
             || seq_results[i].y != mlipp_results[i].y)
                printf("Not equal, %d\n", i);
        }
    }
    else
        printf("Range: n_seq = %d, n_mlipp = %d\n", int(seq_results.size()), result_size);
}

void
knn_test(int n, int k)
{
    MLIPP_KD<int> mlipp;

    vector<Point> data;
    for (int i = 0; i < n; ++i)
        data.push_back((Point){rand(), rand()});

    Point query_point = (Point){rand(), rand()};

    printf("Query point (%d, %d)\n", query_point.x, query_point.y);

    mlipp.bulk_load(data.data(), data.size());

    std::vector<Point> seq_results(k);
    std::vector<double> seq_distances(k);
    int seq_result_size = seq_knn_query(data, query_point, k, 
        seq_results.data(), seq_distances.data());
    std::vector<Point> mlipp_results(k);
    std::vector<double> distances(k);
    int result_size = 0;

    printf("Total size = %d, Result size = %d\n", n, int(seq_results.size()));

    auto start_time = std::chrono::high_resolution_clock::now();

    result_size = mlipp.knn_query(query_point, k, mlipp_results.data(), distances.data());

    auto end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count() * 1e-9;
    cout << "Duration: knn = " << duration << endl;
    cout << "Duration / size: knn = " << duration / result_size << endl;

    printf("result_size = %d = %d\n", seq_result_size, result_size);

    // qsort(seq_results.data(), seq_result_size, sizeof(Point), &compare_pt);
    // qsort(mlipp_results.data(), result_size, sizeof(Point), &compare_pt);

    if (seq_result_size == result_size)
    {
        for (int i = 0; i < seq_result_size; ++i)
        {
            if (seq_results[i].x != mlipp_results[i].x 
             || seq_results[i].y != mlipp_results[i].y) {
                printf("Not equal, %d: seq(%d, %d), mlipp(%d, %d)\n", i, 
                    seq_results[i].x, seq_results[i].y, mlipp_results[i].x, mlipp_results[i].y);
                printf("Distances: seq = %lf, mlipp = %lf\n", seq_distances[i], distances[i]);
            }
        }
    }
    else
        printf("Range: n_seq = %d, n_mlipp = %d\n", seq_result_size, result_size);

    int a;
    double b;
    mlipp.print_depth(&a, &b);
    printf("Tree depth = %d\n", a);
}

int main()
{
    srand(time(NULL));

    run(8*1e7);
    // for (int n = 5e6; n < 1e8; n += 5e6)
    //     run(n);

    // range_test(10000000);

    // knn_test(1000000, 50);

    return 0;
}
