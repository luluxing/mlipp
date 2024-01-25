#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <chrono>

#include "mlipp_zorder.h"
#include "point.h"

using namespace std;

vector<Point> seqRangeQuery(vector<Point> points, 
    Point min_key, Point max_key)
{
    vector<Point> result;
    for (uint i = 0; i < points.size(); ++i)
    {
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
    MLIPP_Z<int> lipp_insert;
    MLIPP_Z<int> lipp_bulk;

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
    MLIPP_Z<int> mlipp_z;

    vector<Point> data;
    for (int i = 0; i < n; ++i)
        data.push_back((Point){rand(), rand()});

    // for (int i = 0; i < n; ++i)
    //     mlipp_z.insert(data[i]);

    mlipp_z.bulk_load(data.data(), data.size());

    int x1 = rand();
    int x2 = rand();
    int y1 = rand();
    int y2 = rand();
    Point min_point = (Point){std::min(x1, x2), std::min(y1, y2)};
    Point max_point = (Point){std::max(x1, x2), std::max(y1, y2)};

    std::vector<Point> seq_result = seqRangeQuery(data, min_point, max_point);
    std::vector<Point> mlipp_z_result(n);
    int result_size = 0;

    printf("Total size = %d, Result size = %d\n", n, int(seq_result.size()));

    auto start_time = std::chrono::high_resolution_clock::now();

    result_size = mlipp_z.range_query(min_point, max_point, mlipp_z_result.data());

    auto end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count() * 1e-9;
    cout << "Duration: range = " << duration << endl;
    cout << "Duration / size: range = " << duration / result_size << endl;


    qsort(seq_result.data(), seq_result.size(), sizeof(Point), &compare_pt);
    qsort(mlipp_z_result.data(), result_size, sizeof(Point), &compare_pt);

    if (int(seq_result.size()) == result_size)
    {
        for (uint i = 0; i < seq_result.size(); ++i)
        {
            if (seq_result[i].x != mlipp_z_result[i].x 
             || seq_result[i].y != mlipp_z_result[i].y)
                printf("Not equal, %d\n", i);
        }
    }
    else
        printf("Range: n_seq = %d, n_mlipp_z = %d\n", int(seq_result.size()), result_size);
}

int main()
{
    srand(time(NULL));

    // run(1e6);
    // for (int n = 5e6; n < 1e8; n += 5e6)
    //     run(n);

    range_test(1000000);

    return 0;
}
