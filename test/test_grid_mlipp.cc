#include <gtest/gtest.h>
#include "grid_mlipp.h"

TEST(GridMlippTest, Construction) {
  GridMlipp<int> mlipp(3, 3, Partition::SPACE);
}

TEST(GridMlippTest, BulkLoadAndExist1) {
  Point<double> points[] = {{1., 2.}, {3., 4.}, {5., 6.}, {7., 8.}, {9., 10.}};
  GridMlipp<double> mlipp(3, 3, Partition::SPACE);
  mlipp.bulk_load(points, 5);
  // Check if the points are correctly partitioned
  EXPECT_TRUE(mlipp.exists({1, 2}));
  EXPECT_TRUE(mlipp.exists({3, 4}));
  EXPECT_FALSE(mlipp.exists({11, 12}));
}

TEST(GridMlippTest, BulkLoadAndExist2) {
  int num = 10000;
  std::vector<Point<double>> data;
  for (int i = 0; i < num; ++i) {
    data.push_back((Point<double>)
          {rand() / static_cast<double>(RAND_MAX),
          rand() / static_cast<double>(RAND_MAX)});
  }
  GridMlipp<double> mlipp(3, 3, Partition::SPACE);
  mlipp.bulk_load(data.data(), data.size());
  // Check if the points are correctly partitioned
  for (const auto& point : data) {
    EXPECT_TRUE(mlipp.exists(point));
  }
  // Check for a point that does not exist
  EXPECT_FALSE(mlipp.exists({1.1, 2.2}));
}

TEST(GridMlippTest, RangeQuery1) {
  Point<double> points[] = {{1., 2.}, {3., 4.}, {5., 6.}, {7., 8.}, {9., 10.}};
  GridMlipp<double> mlipp(3, 3, Partition::SPACE);
  mlipp.bulk_load(points, 5);
  std::vector<Point<double>> result;
  mlipp.rangeQuery({2., 3.}, {8., 9.}, result);
  // Check if the range query returns the correct points
  EXPECT_EQ(result.size(), 3);
  EXPECT_TRUE(std::find(result.begin(), result.end(), Point<double>{3., 4.}) != result.end());
  EXPECT_TRUE(std::find(result.begin(), result.end(), Point<double>{5., 6.}) != result.end());
  EXPECT_TRUE(std::find(result.begin(), result.end(), Point<double>{7., 8.}) != result.end());
  EXPECT_FALSE(std::find(result.begin(), result.end(), Point<double>{9., 10.}) != result.end());
  EXPECT_FALSE(std::find(result.begin(), result.end(), Point<double>{1., 2.}) != result.end());
}

TEST(GridMlippTest, RangeQuery2) {
  int num = 10000;
  int num_queries = 100;
  std::vector<Point<double>> data;
  for (int i = 0; i < num; ++i) {
    data.push_back((Point<double>)
          {rand() / static_cast<double>(RAND_MAX), rand() / static_cast<double>(RAND_MAX)});
  }
  // generate a random range query
  for (int i = 0; i < num_queries; i++) {
    double min_x = rand() / static_cast<double>(RAND_MAX);
    double min_y = rand() / static_cast<double>(RAND_MAX);
    double max_x = min_x + rand() / static_cast<double>(RAND_MAX);
    double max_y = min_y + rand() / static_cast<double>(RAND_MAX);
    Point<double> min_point = {min_x, min_y};
    Point<double> max_point = {max_x, max_y};

    GridMlipp<double> mlipp(3, 3, Partition::SPACE);
    mlipp.bulk_load(data.data(), data.size());
    
    std::vector<Point<double>> result;
    mlipp.rangeQuery(min_point, max_point, result);
    // Check if the range query returns the correct points
    for (const auto& point : data) {
      if (point.x >= min_x && point.x <= max_x &&
          point.y >= min_y && point.y <= max_y) {
        EXPECT_TRUE(std::find(result.begin(), result.end(), point) != result.end());
      } else {
        EXPECT_FALSE(std::find(result.begin(), result.end(), point) != result.end());
      }
    }
  }
}
