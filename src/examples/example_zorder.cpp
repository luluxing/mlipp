#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <chrono>

#include "mlipp_zorder.h"
#include "point.h"

using namespace std;

vector<Point<int>> seq_range_query(vector<Point<int>> points, 
    Point<int> min_key, Point<int> max_key)
{
    vector<Point<int>> result;
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
    MLIPP_Z<int> mlipp;

    /* Generate random data */
    vector<Point<int>> data;
    for (int i = 0; i < n; ++i)
        data.push_back((Point<int>){rand(), rand()});

    /* Generate a random range */
    int x1 = rand();
    int x2 = rand();
    int y1 = rand();
    int y2 = rand();
    Point<int> min_point = (Point<int>){min(x1, x2), min(y1, y2)};
    Point<int> max_point = (Point<int>){max(x1, x2), max(y1, y2)};

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
    vector<Point<int>> seq_result = seq_range_query(data, min_point, max_point);
    vector<Point<int>> mlipp_result(seq_result.size());
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
    run(5e6);
    return 0;
}
