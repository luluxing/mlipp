#include <gtest/gtest.h>
#include <climits>
#include <algorithm>
#include "mbtree.h"

TEST(MBTreeTest, Construction) {
  MBTree<int> mbtree;
  EXPECT_EQ(mbtree.index_size(), 0);
}

TEST(MBTreeTest, BulkLoadSingleLeaf) {
  MBTree<int> mbtree;
  int num = 10;
  std::vector<Point<int>> data;
  for (int i = 0; i < num; ++i) {
    data.push_back((Point<int>){i, i});
  }
  mbtree.bulk_load(data.data(), data.size());
  for (int i = 0; i < num; ++i) {
    EXPECT_TRUE(mbtree.exists(data[i]));
  }

  MBTree<int> mbtree2;
  for (int i = 0; i < num; ++i) {
    mbtree2.insert(data[i]);
  }
  for (int i = 0; i < num; ++i) {
    EXPECT_TRUE(mbtree2.exists(data[i]));
  }
}

TEST(MBTreeTest, BulkLoadTwoLevel) {
  MBTree<int> mbtree;
  int num = 2500;
  std::vector<Point<int>> data;
  for (int i = 0; i < num; ++i) {
    data.push_back((Point<int>){i, i});
  }
  mbtree.bulk_load(data.data(), data.size());
  for (int i = 0; i < num; ++i) {
    EXPECT_TRUE(mbtree.exists(data[i]));
  }

  MBTree<int> mbtree2;
  for (int i = 0; i < num; ++i) {
    mbtree2.insert(data[i]);
  }
  for (int i = 0; i < num; ++i) {
    EXPECT_TRUE(mbtree2.exists(data[i]));
  }
}

TEST(MBTreeTest, BulkLoadRandomData) {
  MBTree<int> mbtree;
  int num = 12500;
  std::vector<Point<int>> data;
  for (int i = 0; i < num; ++i) {
    data.push_back((Point<int>){rand(), rand()});
  }
  mbtree.bulk_load(data.data(), data.size());
  for (int i = 0; i < num; ++i) {
    EXPECT_TRUE(mbtree.exists(data[i]));
  }

  MBTree<int> mbtree2;
  for (int i = 0; i < num; ++i) {
    mbtree2.insert(data[i]);
  }
  for (int i = 0; i < num; ++i) {
    EXPECT_TRUE(mbtree2.exists(data[i]));
  }
}

TEST(MBTreeTest, RangeQueryAll) {
  MBTree<int> mbtree;
  int num = 12500;
  std::vector<Point<int>> data;
  int min_x = INT_MAX, min_y = INT_MAX, max_x = INT_MIN, max_y = INT_MIN;
  for (int i = 0; i < num; ++i) {
    data.push_back((Point<int>){rand(), rand()});
    min_x = std::min(min_x, data[i].x);
    min_y = std::min(min_y, data[i].y);
    max_x = std::max(max_x, data[i].x);
    max_y = std::max(max_y, data[i].y);
  }
  // Remove the duplicate points
  data.erase(std::unique(data.begin(), data.end()), data.end());
  num = data.size();

  mbtree.bulk_load(data.data(), data.size());
  std::vector<Point<int>> result;
  mbtree.rangeQuery((Point<int>){min_x, min_y}, (Point<int>){max_x, max_y}, result);
  EXPECT_EQ(result.size(), num);
  // Sort both data and result and check if they are the same
  std::sort(data.begin(), data.end());
  std::sort(result.begin(), result.end());
  for (int i = 0; i < num; ++i) {
    EXPECT_EQ(data[i], result[i]);
  }
  
  MBTree<int> mbtree2;
  for (int i = 0; i < num; ++i) {
    mbtree2.insert(data[i]);
  }
  std::vector<Point<int>> result2;
  mbtree2.rangeQuery((Point<int>){min_x, min_y}, (Point<int>){max_x, max_y}, result2);
  EXPECT_EQ(result2.size(), num);
  // Sort both data and result and check if they are the same
  std::sort(data.begin(), data.end());
  std::sort(result2.begin(), result2.end());
  for (int i = 0; i < num; ++i) {
    EXPECT_EQ(data[i], result2[i]);
  }
}

TEST(MBTreeTest, RangeQueryPartial) {
  MBTree<int> mbtree;
  int num = 12500;
  int query_num = 1000;
  std::vector<Point<int>> data;
  int min_x = INT_MAX, min_y = INT_MAX, max_x = INT_MIN, max_y = INT_MIN;
  for (int i = 0; i < num; ++i) {
    data.push_back((Point<int>){rand(), rand()});
    min_x = std::min(min_x, data[i].x);
    min_y = std::min(min_y, data[i].y);
    max_x = std::max(max_x, data[i].x);
    max_y = std::max(max_y, data[i].y);
  }
  // Remove the duplicate points
  data.erase(std::unique(data.begin(), data.end()), data.end());
  num = data.size();
  // Generate rectangle query
  std::vector<std::pair<Point<int>, Point<int>>> query_data;
  for (int i = 0; i < query_num; ++i) {
    // Randomly pick a point between min_x and max_x, min_y and max_y
    int x1 = min_x + rand() % (max_x - min_x + 1);
    int y1 = min_y + rand() % (max_y - min_y + 1);
    int x2 = min_x + rand() % (max_x - min_x + 1);
    int y2 = min_y + rand() % (max_y - min_y + 1);
    // Find the lower left point and upper right point
    int lower_x = std::min(x1, x2);
    int lower_y = std::min(y1, y2);
    int upper_x = std::max(x1, x2);
    int upper_y = std::max(y1, y2);
    query_data.push_back(std::make_pair((Point<int>){lower_x, lower_y}, (Point<int>){upper_x, upper_y}));
  }

  // mbtree.bulk_load(data.data(), data.size());
  for (int i = 0; i < num; ++i) {
    mbtree.insert(data[i]);
  }

  for (int i = 0; i < query_num; ++i) {
    std::vector<Point<int>> result;
    // Range query
    mbtree.rangeQuery(query_data[i].first, query_data[i].second, result);
    // Linearly scan the data and find all the matching points
    std::vector<Point<int>> result2;
    for (int j = 0; j < num; ++j) {
      if (data[j].x >= query_data[i].first.x && data[j].x <= query_data[i].second.x &&
          data[j].y >= query_data[i].first.y && data[j].y <= query_data[i].second.y) {
        result2.push_back(data[j]);
      }
    }
    // Sort both result and result2 and check if they are the same
    std::sort(result.begin(), result.end());
    std::sort(result2.begin(), result2.end());
    EXPECT_EQ(result.size(), result2.size());
    for (int j = 0; j < result.size(); ++j) {
      EXPECT_EQ(result[j], result2[j]);
    }
  }
}