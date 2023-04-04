#include <lipp.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <chrono>
#include <algorithm>

using namespace std;

typedef pair<int, int> Point;

bool sortbyfirst(const Point &a,
              const Point &b)
{
    return (a.first < b.first);
}

vector<Point> seqRangeQuery(vector<Point> points, 
    Point min_key, Point max_key)
{
    vector<Point> result;
    for (uint i = 0; i < points.size(); ++i)
    {
        if (points[i].first >= min_key.first 
            && points[i].first <= max_key.first)
            result.push_back(points[i]);
    }
    return result;
}

void
run(int n_loop)
{
    LIPP<int, int> lipp_insert;
    LIPP<int, int> lipp_bulk;

    vector<Point> data_raw;
    for (int i = 0; i < n_loop; ++i)
        data_raw.push_back(make_pair(rand(), rand()));

    sort(data_raw.begin(), data_raw.end(), &sortbyfirst);

    vector<int> idx;
    idx.push_back(0);
    for (int i = 1; i < n_loop; ++i)
        if (data_raw[i].first > data_raw[i - 1].first)
            idx.push_back(i);

    vector<Point> data;
    for (int i = 0; i < int(idx.size()); ++i)
        data.push_back(data_raw[idx[i]]);

    random_shuffle(data.begin(), data.end());

    int n = data.size();

    auto start_time = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < n; ++i)
        lipp_insert.insert(data[i]);

    auto end_time = chrono::high_resolution_clock::now();
    auto duration_insert = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count() * 1e-9;

    start_time = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < n; ++i)
        lipp_insert.exists(data[i].first);

    end_time = chrono::high_resolution_clock::now();
    auto duration_scan_insert = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count() * 1e-9;

    sort(data.begin(), data.end(), &sortbyfirst);

    start_time = std::chrono::high_resolution_clock::now();

    lipp_bulk.bulk_load(data.data(), data.size());

    end_time = chrono::high_resolution_clock::now();
    auto duration_build = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count() * 1e-9;

    start_time = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < n; ++i)
        lipp_bulk.exists(data[i].first);

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
range_test(int n_loop)
{
    LIPP<int, int> lipp;

    vector<Point> data_raw;
    for (int i = 0; i < n_loop; ++i)
        data_raw.push_back(make_pair(rand(), rand()));

    sort(data_raw.begin(), data_raw.end(), &sortbyfirst);

    vector<int> idx;
    idx.push_back(0);
    for (int i = 1; i < n_loop; ++i)
        if (data_raw[i].first > data_raw[i - 1].first)
            idx.push_back(i);

    vector<Point> data;
    for (int i = 0; i < int(idx.size()); ++i)
        data.push_back(data_raw[idx[i]]);

    lipp.bulk_load(data.data(), data.size());

    Point xmin = make_pair(rand(), rand());
    Point xmax;
    do {
        xmax = make_pair(rand(), rand());
    } while (xmax.first <= xmin.first);

    std::vector<Point> seq_result = seqRangeQuery(data, xmin, xmax);
    std::vector<Point> lipp_result = lipp.rangeQuery(xmin.first, xmax.first);

    printf("n_seq = %d, n_lipp = %d\n", int(seq_result.size()), int(lipp_result.size()));
    if (seq_result.size() == lipp_result.size())
    {
        for (uint i = 0; i < seq_result.size(); ++i)
        {
            if (seq_result[i].first != lipp_result[i].first)
                printf("Not equal, %d\n", i);
        }
    }
}

int main()
{
    srand(time(NULL));

    run(1e6);
    for (int n = 5e6; n < 1e8; n += 5e6)
        run(n);

    // range_test(1000000);

    return 0;
}
