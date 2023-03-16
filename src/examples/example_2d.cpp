#include <lipp_2d.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <chrono>

using namespace std;

typedef pair<int, int> Point;

int main()
{
    srand(time(NULL));

    LIPP<int> lipp;

    int n = 10000000;
    Point p1 = make_pair(104394, 2382034);

    auto start_time = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < n / 2; ++i)
    {
        Point p = make_pair(rand(), rand());
        // if (!lipp.exists(p))
        lipp.insert(p);
    }

    if (!lipp.exists(p1))
        lipp.insert(p1);
    else
        cout << "Point already exists!" << endl;

    for (int i = 0; i < n / 2; ++i)
    {
        Point p = make_pair(rand(), rand());
        // if (!lipp.exists(p))
        lipp.insert(p);
    }

    cout << "Point found = " << (lipp.exists(p1) ? "true" : "false") << endl;

    auto end_time = chrono::high_resolution_clock::now();
    auto duration = end_time - start_time;
    auto duration_seconds = chrono::duration_cast<chrono::nanoseconds>(duration).count() * 1e-9;

    cout << "Duration = " << duration_seconds << "s" << endl;
    cout << "Throughput = " << n / (1000000 * duration_seconds) << "M/s" << endl;

    // lipp.verify();
    // lipp.print_depth();
    // lipp.print_stats();

    return 0;
}
