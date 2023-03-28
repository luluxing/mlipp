#include <lipp_2darray.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <chrono>

using namespace std;

typedef pair<int, int> Point;

int main()
{
    srand(time(NULL));

    for (int n = 1e6; n < 1e8; n += 1e6)
    {
        LIPP<int> lipp_insert;
        LIPP<int> lipp_bulk;

        vector<Point> data;
        for (int i = 0; i < n; ++i)
            data.push_back(make_pair(rand(), rand()));

        auto start_time = std::chrono::high_resolution_clock::now();

        /*for (int i = 0; i < n; ++i)
            lipp_insert.insert(data[i]);*/

        auto end_time = chrono::high_resolution_clock::now();
        auto duration_insert = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count() * 1e-9;

        start_time = std::chrono::high_resolution_clock::now();

        /*for (int i = 0; i < n; ++i)
            lipp_insert.exists(data[i]);*/

        end_time = chrono::high_resolution_clock::now();
        auto duration_scan_insert = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count() * 1e-9;

        start_time = std::chrono::high_resolution_clock::now();

        lipp_bulk.bulk_load(data.data(), data.size());

        end_time = chrono::high_resolution_clock::now();
        auto duration_build = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count() * 1e-9;

        start_time = std::chrono::high_resolution_clock::now();

        for (int i = 0; i < n; ++i)
            lipp_bulk.exists(data[i]);

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

    return 0;
}
