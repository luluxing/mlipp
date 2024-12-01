#include <gtest/gtest.h>
#include <climits>
#include "mlipp_kd.h"
#include "point.h"

int num = 100000;

TEST(MLIPPTest, InitialConstructionAndSingleInsertion) {
  MLIPP_KD<int> mlipp;
  mlipp.insert((Point){10, 20});

  Point min_point = (Point){5, 10};
  Point max_point = (Point){15, 30};

  std::vector<Point> result;
  mlipp.rangeQuery(min_point, max_point, result);
  EXPECT_EQ(result.size(), 1);
  EXPECT_EQ(result[0].x, 10);
  EXPECT_EQ(result[0].y, 20);
}

TEST(MLIPPTest, BulkloadingAndExistingData) {
  MLIPP_KD<int> mlipp;
  std::vector<Point> data;
  for (int i = 0; i < num; ++i) {
    data.push_back((Point){rand(), rand()});
  }
  mlipp.bulk_load(data.data(), data.size());
  // Test existence
  for (int i = 0; i < num; ++i) {
    EXPECT_TRUE(mlipp.exists(data[i]));
  }
}

TEST(MLIPPTest, NonExistingData) {
  MLIPP_KD<int> mlipp;
  std::vector<Point> data;
  for (int i = 0; i < num; ++i) {
    data.push_back((Point){rand(), rand()});
  }
  mlipp.bulk_load(data.data(), data.size());
  // Test non-existence
  for (int i = 0; i < num; ++i) {
    int x = rand();
    int y = rand();
    if (std::find(data.begin(), data.end(), (Point){x, y}) == data.end()) {
      EXPECT_FALSE(mlipp.exists((Point){x, y}));
    }
  }
}

TEST(MLIPPTest, RangeQuery) {
  MLIPP_KD<int> mlipp;
  std::vector<Point> data;
  int min_x = INT_MAX, min_y = INT_MAX, max_x = INT_MIN, max_y = INT_MIN;
  for (int i = 0; i < num; ++i) {
    data.push_back((Point){rand(), rand()});
    min_x = std::min(min_x, data[i].x);
    min_y = std::min(min_y, data[i].y);
    max_x = std::max(max_x, data[i].x);
    max_y = std::max(max_y, data[i].y);
  }
  // Bulk load
  mlipp.bulk_load(data.data(), data.size());
  // Range query
  std::vector<Point> result;
  mlipp.rangeQuery((Point){min_x, min_y}, (Point){max_x, max_y}, result);
  EXPECT_EQ(result.size(), num);
  for (int i = 0; i < num; ++i) {
    EXPECT_TRUE(std::find(result.begin(), result.end(), data[i]) != result.end());
  }
}