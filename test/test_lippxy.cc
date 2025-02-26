#include <gtest/gtest.h>
#include <climits>
#include "two_lipp.h"
#include "point.h"

TEST(LIPPXYTest, BulkloadAndExist) {
  {
    LIPP_XY<int> mlipp;
    int num = 100;
    std::vector<Point<int>> data;
    for (int i = 0; i < num; ++i) {
      data.push_back((Point<int>){i * 3, i * 4});
    }
    mlipp.bulk_load(data.data(), data.size());
    // Test existence
    for (int i = 0; i < num; ++i) {
      EXPECT_TRUE(mlipp.exists(data[i]));
    }
    // Test non-existence
    for (int i = 0; i < num; ++i) {
      EXPECT_FALSE(mlipp.exists((Point<int>){i * 3 + 1, i * 4 + 1}));
    }
  }
  {
    LIPP_XY<double> mlipp;
    int num = 100;
    std::vector<Point<double>> data;
    for (int i = 0; i < num; ++i) {
      data.push_back((Point<double>)
        {i * 3.0 / static_cast<double>(RAND_MAX), i * 4.0 / static_cast<double>(RAND_MAX)});
    }
    mlipp.bulk_load(data.data(), data.size());
    // Test existence
    for (int i = 0; i < num; ++i) {
      EXPECT_TRUE(mlipp.exists(data[i]));
    }
    // Test non-existence
    for (int i = 0; i < num; ++i) {
      EXPECT_FALSE(mlipp.exists((Point<double>)
        {i * 3.0 / static_cast<double>(RAND_MAX) + 1, i * 4.0 / static_cast<double>(RAND_MAX) + 1}));
    }
  }
}

TEST(LIPPXYTest, InitialConstructionAndSingleInsertion) {
  
}

TEST(LIPPXYTest, BulkloadingAndExistingData) {
  
}

TEST(LIPPXYTest, NonExistingData) {
  
}

TEST(LIPPXYTest, RangeQuery) {
  
}

TEST(LIPPXYTest, DuplicateCoordinates) {

}