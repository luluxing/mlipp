#include <gtest/gtest.h>
#include "point.h"

TEST(PointTest, PointEquality) {
  Point<int> p1 = {1, 2};
  Point<int> p2 = {1, 2};
  Point<int> p3 = {2, 3};
  EXPECT_TRUE(p1 == p2);
  EXPECT_FALSE(p1 == p3);
}

TEST(PointTest, PointComparison) {
  Point<int> p1 = {1, 5};
  Point<int> p2 = {1, 5};
  Point<int> p3 = {2, 3};
  EXPECT_EQ(compare_x<int>(&p1, &p2), 0);
  EXPECT_EQ(compare_y<int>(&p1, &p2), 0);
  EXPECT_EQ(compare_pt<int>(&p1, &p2), 0);
  EXPECT_EQ(compare_x<int>(&p1, &p3), -1);
  EXPECT_EQ(compare_y<int>(&p1, &p3), 1);
  EXPECT_EQ(compare_pt<int>(&p1, &p3), -1);
}

TEST(PointTest, PointValue) {
  Point<int> p = {1, 2};
  EXPECT_EQ(PT_VAL(p, 0), 1);
  EXPECT_EQ(PT_VAL(p, 1), 2);
  EXPECT_TRUE(PT_EQ(p, p));

  Point<double> p2 = {1.0, 2.0};
  EXPECT_EQ(PT_VAL(p2, 0), 1.0);
  EXPECT_EQ(PT_VAL(p2, 1), 2.0);
  EXPECT_TRUE(PT_EQ(p2, p2));
}