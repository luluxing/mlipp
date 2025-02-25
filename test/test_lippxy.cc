#include <gtest/gtest.h>
#include <climits>
#include "two_lipp.h"
#include "point.h"

int num = 100000;

TEST(LIPPXYTest, InitialConstructionAndSingleInsertion) {
  {
    LIPP_XY<int> mlipp;
    mlipp.insert((Point<int>){10, 20});

    Point<int> min_point = (Point<int>){5, 10};
    Point<int> max_point = (Point<int>){15, 30};

    std::vector<Point<int>> result;
    mlipp.rangeQuery(min_point, max_point, result);
    EXPECT_EQ(result.size(), 1);
    EXPECT_EQ(result[0].x, 10);
    EXPECT_EQ(result[0].y, 20);
  }
  {
    LIPP_XY<double> mlipp;
    mlipp.insert((Point<double>){10.5, 20.5});

    Point<double> min_point = (Point<double>){5.5, 10.5};
    Point<double> max_point = (Point<double>){15.5, 30.5};

    std::vector<Point<double>> result;
    mlipp.rangeQuery(min_point, max_point, result);
    EXPECT_EQ(result.size(), 1);
    EXPECT_EQ(result[0].x, 10.5);
    EXPECT_EQ(result[0].y, 20.5);
  }
}

TEST(LIPPXYTest, BulkloadingAndExistingData) {
  {
    LIPP_XY<int> mlipp;
    std::vector<Point<int>> data;
    for (int i = 0; i < num; ++i) {
      data.push_back((Point<int>){rand(), rand()});
    }
    mlipp.bulk_load(data.data(), data.size());
    // Test existence
    for (int i = 0; i < num; ++i) {
      EXPECT_TRUE(mlipp.exists(data[i]));
    }
  }
  {
    LIPP_XY<double> mlipp;
    std::vector<Point<double>> data;
    for (int i = 0; i < num; ++i) {
      data.push_back((Point<double>)
        {rand() / static_cast<double>(RAND_MAX), rand() / static_cast<double>(RAND_MAX)});
    }
    mlipp.bulk_load(data.data(), data.size());
    // Test existence
    for (int i = 0; i < num; ++i) {
      EXPECT_TRUE(mlipp.exists(data[i]));
    }
  }
}

TEST(LIPPXYTest, NonExistingData) {
  int num = 10000;
  {
    LIPP_XY<int> mlipp;
    std::vector<Point<int>> data;
    for (int i = 0; i < num; ++i) {
      data.push_back((Point<int>){rand(), rand()});
    }
    mlipp.bulk_load(data.data(), data.size());
    // Test non-existence
    for (int i = 0; i < num / 10; ++i) {
      int x = rand();
      int y = rand();
      if (std::find(data.begin(), data.end(), (Point<int>){x, y}) == data.end()) {
        EXPECT_FALSE(mlipp.exists((Point<int>){x, y}));
      }
    }
  }

  {
    LIPP_XY<double> mlipp;
    std::vector<Point<double>> data;
    for (int i = 0; i < num; ++i) {
      data.push_back((Point<double>)
        {rand() / static_cast<double>(RAND_MAX), rand() / static_cast<double>(RAND_MAX)});
    }
    mlipp.bulk_load(data.data(), data.size());
    // Test non-existence
    for (int i = 0; i < num / 10; ++i) {
      double x = rand() / static_cast<double>(RAND_MAX);
      double y = rand() / static_cast<double>(RAND_MAX);
      if (std::find(data.begin(), data.end(), (Point<double>){x, y}) == data.end()) {
        EXPECT_FALSE(mlipp.exists((Point<double>){x, y}));
      }
    }
  }
}

TEST(LIPPXYTest, RangeQuery) {
  {
    LIPP_XY<int> mlipp;
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
    mlipp.bulk_load(data.data(), data.size());
    // Range query
    std::vector<Point<int>> result;
    mlipp.rangeQuery((Point<int>){min_x, min_y}, (Point<int>){max_x, max_y}, result);
    EXPECT_EQ(result.size(), num);
    // Sort both data and result and check if they are the same
    std::sort(data.begin(), data.end());
    std::sort(result.begin(), result.end());
    for (int i = 0; i < num; ++i) {
      EXPECT_EQ(data[i], result[i]);
    }
  }

  {
    LIPP_XY<double> mlipp;
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
    mlipp.bulk_load(data.data(), data.size());
    // Range query
    std::vector<Point<double>> result;
    mlipp.rangeQuery((Point<double>){min_x, min_y}, (Point<double>){max_x, max_y}, result);
    EXPECT_EQ(result.size(), num);
    // Sort both data and result and check if they are the same
    std::sort(data.begin(), data.end());
    std::sort(result.begin(), result.end());
    for (int i = 0; i < num; ++i) {
      EXPECT_EQ(data[i], result[i]);
    }
  }
}