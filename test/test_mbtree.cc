#include <gtest/gtest.h>
#include <climits>
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
}