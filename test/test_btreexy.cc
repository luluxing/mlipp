#include <gtest/gtest.h>
#include <climits>
#include "btree_xy.h"
#include "point.h"

TEST(BTreeXYTest, BulkloadAndExist) {
  {
    BTreeXY<int> btreexy;
    int num = 10000;
    std::vector<Point<int>> data;
    for (int i = 0; i < num; ++i) {
      data.push_back((Point<int>){i * 3, i * 4});
    }
    btreexy.bulk_load(data.data(), data.size());
    // Test existence
    for (int i = 0; i < num; ++i) {
      EXPECT_TRUE(btreexy.exists(data[i]));
    }
    // Test non-existence
    for (int i = 0; i < num; ++i) {
      EXPECT_FALSE(btreexy.exists((Point<int>){i * 3 + 1, i * 4 + 1}));
    }
  }
  {
    BTreeXY<double> btreexy;
    int num = 100;
    std::vector<Point<double>> data;
    for (int i = 0; i < num; ++i) {
      data.push_back((Point<double>)
        {i * 3.0 / static_cast<double>(RAND_MAX), i * 4.0 / static_cast<double>(RAND_MAX)});
    }
    btreexy.bulk_load(data.data(), data.size());
    // Test existence
    for (int i = 0; i < num; ++i) {
      EXPECT_TRUE(btreexy.exists(data[i]));
    }
    // Test non-existence
    for (int i = 0; i < num; ++i) {
      EXPECT_FALSE(btreexy.exists((Point<double>)
        {i * 3.0 / static_cast<double>(RAND_MAX) + 1, i * 4.0 / static_cast<double>(RAND_MAX) + 1}));
    }
  }
}

TEST(BTreeXYTest, RangeQuery) {
  int num = 100000;
  {
    BTreeXY<int> btreexy;
    std::vector<Point<int>> data;
    int min_x = INT_MAX, min_y = INT_MAX, max_x = INT_MIN, max_y = INT_MIN;
    for (int i = 0; i < num; ++i) {
      data.push_back((Point<int>){rand(), rand()});
      min_x = std::min(min_x, data[i].x);
      min_y = std::min(min_y, data[i].y);
      max_x = std::max(max_x, data[i].x);
      max_y = std::max(max_y, data[i].y);
    }
    // Bulk load
    btreexy.bulk_load(data.data(), data.size());
    // Range query
    std::vector<Point<int>> result;
    btreexy.rangeQuery((Point<int>){min_x, min_y}, (Point<int>){max_x, max_y}, result);
    EXPECT_EQ(result.size(), num);
    // Sort both data and result and check if they are the same
    std::sort(data.begin(), data.end());
    std::sort(result.begin(), result.end());
    for (int i = 0; i < num; ++i) {
      EXPECT_EQ(data[i], result[i]);
    }
  }
  {
    BTreeXY<double> btreexy;
    std::vector<Point<double>> data;
    double min_x = 1.0, min_y = 1.0, max_x = 0.0, max_y = 0.0;
    for (int i = 0; i < num; ++i) {
      data.push_back((Point<double>)
        {rand() / static_cast<double>(RAND_MAX), rand() / static_cast<double>(RAND_MAX)});
      min_x = std::min(min_x, data[i].x);
      min_y = std::min(min_y, data[i].y);
      max_x = std::max(max_x, data[i].x);
      max_y = std::max(max_y, data[i].y);
    }
    // Bulk load
    btreexy.bulk_load(data.data(), data.size());
    // Range query
    std::vector<Point<double>> result;
    btreexy.rangeQuery((Point<double>){min_x, min_y}, (Point<double>){max_x, max_y}, result);
    EXPECT_EQ(result.size(), num);
    // Sort both data and result and check if they are the same
    std::sort(data.begin(), data.end());
    std::sort(result.begin(), result.end());
    for (int i = 0; i < num; ++i) {
      EXPECT_EQ(data[i], result[i]);
    }
  }
}

TEST(BTreeXYTest, DuplicateCoordinates) {
  int num = 10;
  {
    BTreeXY<int> btreexy;
    std::vector<Point<int>> data;
    int min_x = INT_MAX, min_y = INT_MAX, max_x = INT_MIN, max_y = INT_MIN;
    for (int i = 0; i < num; ++i) {
      data.push_back((Point<int>){i+1, i+1});
    }
    data.push_back((Point<int>){1, 2});
    data.push_back((Point<int>){1, 3});

    // Bulk load
    btreexy.bulk_load(data.data(), data.size());

    // Existence
    EXPECT_TRUE(btreexy.exists((Point<int>){1, 2}));
    EXPECT_TRUE(btreexy.exists((Point<int>){1, 3}));
    EXPECT_FALSE(btreexy.exists((Point<int>){1, 4}));

    // Range query: skinny rectange
    std::vector<Point<int>> result;
    btreexy.rangeQuery((Point<int>){0, 0}, (Point<int>){2, 3}, result);
    EXPECT_EQ(result.size(), 4);
    // Sort both data and result and check if they are the same
    std::sort(data.begin(), data.end());
    std::sort(result.begin(), result.end());
    for (int i = 0; i < result.size(); ++i) {
      EXPECT_EQ(data[i], result[i]);
    }

    // Range query: wide rectange
    result.clear();
    btreexy.rangeQuery((Point<int>){0, 0}, (Point<int>){10, 3}, result);
    EXPECT_EQ(result.size(), 5);
    // Sort both data and result and check if they are the same
    std::sort(data.begin(), data.end());
    std::sort(result.begin(), result.end());
    for (int i = 0; i < result.size(); ++i) {
      EXPECT_EQ(data[i], result[i]);
    }
  }
}