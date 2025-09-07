#include <gtest/gtest.h>
#include <climits>
#include "mbtree.h"

TEST(MBTreeTest, Construction) {
  MBTree<int> mbtree;
  EXPECT_EQ(mbtree.index_size(), 0);
}