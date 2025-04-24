#include <gtest/gtest.h>
#include "multi_mlipp.h"

TEST(MultiMlippTest, Construction) {
  {
    MultiMlipp<int> mlipp(3, Axis::X_AXIS, Partition::DATA);
  }
  {
    MultiMlipp<int> mlipp(3, Axis::Y_AXIS, Partition::SPACE);
  }
}

TEST(MultiMlippTest, BulkLoadAndExist1) {
  std::vector<Axis> axes = {Axis::X_AXIS, Axis::Y_AXIS};
  std::vector<Partition> partitions = {Partition::DATA, Partition::SPACE};
  for (const auto& axis : axes) {
    for (const auto& partition : partitions) {
      Point<double> points[] = {{1., 2.}, {3., 4.}, {5., 6.}, {7., 8.}, {9., 10.}};
      MultiMlipp<double> mlipp(3, axis, partition);
      mlipp.bulk_load(points, 5);
      // Check if the points are correctly partitioned
      EXPECT_TRUE(mlipp.exists({1, 2}));
      EXPECT_TRUE(mlipp.exists({3, 4}));
      EXPECT_FALSE(mlipp.exists({11, 12}));
    }
  }
}

TEST(MultiMlippTest, BulkLoadAndExist2) {
  int num = 10000;
  std::vector<Axis> axes = {Axis::X_AXIS, Axis::Y_AXIS};
  std::vector<Partition> partitions = {Partition::DATA, Partition::SPACE};
  for (const auto& axis : axes) {
    for (const auto& partition : partitions) {
      std::vector<Point<double>> data;
      for (int i = 0; i < num; ++i) {
        data.push_back((Point<double>)
          {rand() / static_cast<double>(RAND_MAX), rand() / static_cast<double>(RAND_MAX)});
      }
      MultiMlipp<double> mlipp(3, axis, partition);
      mlipp.bulk_load(data.data(), data.size());
      // Check if the points are correctly partitioned
      for (const auto& point : data) {
        EXPECT_TRUE(mlipp.exists(point));
      }
      // Check for a point that does not exist
      // EXPECT_FALSE(mlipp.exists({1.1, 2.2}));
    }
  }
}