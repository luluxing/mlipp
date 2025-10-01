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
  MBTree<uint64_t> mbtree;
  int num = 12500;
  std::vector<Point<uint64_t>> data;
  for (int i = 0; i < num; ++i) {
    data.push_back((Point<uint64_t>){rand(), rand()});
  }
  mbtree.bulk_load(data.data(), data.size());
  for (int i = 0; i < num; ++i) {
    EXPECT_TRUE(mbtree.exists(data[i]));
  }

  MBTree<uint64_t> mbtree2;
  for (int i = 0; i < num; ++i) {
    mbtree2.insert(data[i]);
  }
  for (int i = 0; i < num; ++i) {
    EXPECT_TRUE(mbtree2.exists(data[i]));
  }
}