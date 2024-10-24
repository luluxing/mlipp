#include <iostream>
#include <fstream>
#include <cstring>
#include <stdlib.h>
#include <time.h>
#include <chrono>

#include "lipp.h"

using namespace std;

typedef pair<int, int> Point;

bool sortbyfirst(const Point &a, const Point &b)
{
    return (a.first < b.first);
}

double run_lipp(vector<Point> data)
{
    LIPP<int, int> index;
    printf("Building index with %lu points\n", data.size());
    index.bulk_load(data.data(), data.size());

    printf("Querying index\n");
    auto start_time = chrono::high_resolution_clock::now();
    for (uint i = 0; i < data.size(); ++i)
    {
        if (!index.exists(data[i].first))
            printf("Point not found: %f, %f\n", data[i].first, data[i].second);
    }
    auto end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count() * 1e-9;
    return duration / data.size();
}

vector<Point> get_points(string filename)
{
    vector<Point> result;

    fstream fin;
    fin.open(filename, ios::in);

    // if (fin.fail())
    //     printf("Failed\n");

    // Get the roll number
    // of which the data is required
    int count = 0;
    int lat, lon; 
  
    // Read the Data from the file
    // as String Vector
    vector<string> row;
    string line, word, temp;
  
    while (!fin.eof()) {
  
        row.clear();
        getline(fin, line, '\n');
  
        // used for breaking words
        stringstream s(line);
  
        // read every column data of a row and
        // store it in a string variable, 'word'
        while (getline(s, word, ',')) {
            row.push_back(word);
        }

        lat = stoi(row[0]);
        lon = stoi(row[1]);

        result.push_back(make_pair(lat, lon));
        count += 1;
    }

    printf("Read %d rows\n", count);

    return result;
}

void run(string filename)
{
    vector<Point> data = get_points(filename);
    sort(data.begin(), data.end(), &sortbyfirst);
    printf("Running mlipp_kd\n");
    auto duration = run_lipp(data);
    printf("Time per query: %f\n", duration);
}

int main()
{
    run("/home/maxime/Documents/phd/data/ucr_star/tweets/tweets_data_1d.csv");
    return 0;
}