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