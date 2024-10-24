#include <iostream>
#include <fstream>
#include <cstring>
#include <stdlib.h>
#include <time.h>
#include <chrono>

#include "point.h"
#include "mlipp_kd.h"
#include "mlipp_zorder.h"

using namespace std;

double run_mlipp_kd(vector<Point> data)
{
    MLIPP_KD<float> index;
    printf("Building index with %d points\n", data.size());
    index.bulk_load(data.data(), data.size());

    printf("Querying index\n");
    auto start_time = chrono::high_resolution_clock::now();
    for (uint i = 0; i < data.size(); ++i)
    {
        if (!index.exists(data[i]))
            printf("Point not found: %f, %f\n", data[i].x, data[i].y);
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
    float lat, lon; 
  
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

        lat = stof(row[0]);
        lon = stof(row[1]);

        result.push_back((Point){lat, lon});
        count++;
    }

    printf("Read %d rows\n", count);

    return result;
}

void run(string filename)
{
    vector<Point> data = get_points(filename);
    printf("Running mlipp_kd\n");
    auto duration = run_mlipp_kd(data);
    printf("Time per query: %f\n", duration);
}

int main()
{
    run("/home/maxime/Documents/phd/data/ucr_star/tweets/tweets_data.csv");
    return 0;
}