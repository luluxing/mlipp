#include <gtest/gtest.h>
#include <algorithm>
#include <cstdlib>
#include <vector>

#include "btree_zorder.h"
#include "point.h"

TEST(BTreeZOrderTest, BulkLoadAndLookup) {
  BTreeZOrder<int> tree;
  std::vector<Point<int>> data;
  for (int i = 0; i < 500; ++i) {
    data.push_back((Point<int>){i * 2, i * 3});
  }

  tree.bulk_load(data.data(), static_cast<int>(data.size()));

  for (const auto& pt : data) {
    Point<int> found;
    EXPECT_TRUE(tree.lookup(pt, found));
    EXPECT_EQ(found, pt);
  }

  EXPECT_FALSE(tree.exists((Point<int>){-1, -1}));
}

TEST(BTreeZOrderTest, InsertAndDuplicateKeys) {
  BTreeZOrder<int> tree;
  std::vector<Point<int>> data = {
      {1, 1}, {2, 2}, {3, 3}, {1, 1}, {2, 2}, {5, 7}};

  for (const auto& pt : data) {
    tree.insert(pt);
  }

  for (const auto& pt : data) {
    Point<int> found;
    EXPECT_TRUE(tree.lookup(pt, found));
    EXPECT_EQ(found, pt);
  }

  std::vector<Point<int>> range_result;
  tree.rangeQuery((Point<int>){0, 0}, (Point<int>){10, 10}, range_result);
  ASSERT_EQ(range_result.size(), data.size());

  auto sorted_expected = data;
  std::sort(sorted_expected.begin(), sorted_expected.end());
  std::sort(range_result.begin(), range_result.end());
  EXPECT_EQ(sorted_expected, range_result);
}

TEST(BTreeZOrderTest, RangeQueryMatchesBruteForce) {
  BTreeZOrder<int> tree;
  std::vector<Point<int>> data;
  const int num = 2000;
  srand(0);
  for (int i = 0; i < num; ++i) {
    data.push_back((Point<int>){rand() % 1000, rand() % 1000});
  }
  tree.bulk_load(data.data(), static_cast<int>(data.size()));

  for (int q = 0; q < 5; ++q) {
    int x1 = rand() % 1000;
    int x2 = rand() % 1000;
    int y1 = rand() % 1000;
    int y2 = rand() % 1000;
    Point<int> min_key = {std::min(x1, x2), std::min(y1, y2)};
    Point<int> max_key = {std::max(x1, x2), std::max(y1, y2)};

    std::vector<Point<int>> expected;
    for (const auto& pt : data) {
      if (pt.x >= min_key.x && pt.x <= max_key.x && pt.y >= min_key.y &&
          pt.y <= max_key.y) {
        expected.push_back(pt);
      }
    }

    std::vector<Point<int>> result;
    tree.rangeQuery(min_key, max_key, result);

    std::sort(expected.begin(), expected.end());
    std::sort(result.begin(), result.end());
    EXPECT_EQ(expected, result);
  }
}

TEST(BTreeZOrderTest, SupportsDoubleCoordinates) {
  BTreeZOrder<double> tree;
  std::vector<Point<double>> data;
  for (int i = 0; i < 100; ++i) {
    data.push_back(
        (Point<double>){i * 0.01, (100 - i) * 0.01});
  }
  tree.bulk_load(data.data(), static_cast<int>(data.size()));

  for (const auto& pt : data) {
    Point<double> found;
    EXPECT_TRUE(tree.lookup(pt, found));
    EXPECT_DOUBLE_EQ(found.x, pt.x);
    EXPECT_DOUBLE_EQ(found.y, pt.y);
  }

  std::vector<Point<double>> result;
  tree.rangeQuery((Point<double>){0.2, 0.2}, (Point<double>){0.8, 0.8},
                  result);
  for (const auto& pt : result) {
    EXPECT_GE(pt.x, 0.2);
    EXPECT_LE(pt.x, 0.8);
    EXPECT_GE(pt.y, 0.2);
    EXPECT_LE(pt.y, 0.8);
  }
}
