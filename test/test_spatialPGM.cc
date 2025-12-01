#include <algorithm>
#include <climits>
#include <cstdlib>
#include <vector>

#include <gtest/gtest.h>

#include "spatialPGM.h"

namespace {

template <class T>
void sort_points(std::vector<Point<T>>& points) {
  std::sort(points.begin(), points.end(), [](const Point<T>& lhs, const Point<T>& rhs) {
    if (lhs.x == rhs.x) return lhs.y < rhs.y;
    return lhs.x < rhs.x;
  });
}

}  // namespace

TEST(SpatialPGMTest, BuildsHierarchicalSegments) {
  std::vector<Point<int>> points = {
      {0, 0}, {1, 0}, {2, 0}, {10, 1}, {11, 2}, {12, 2}};

  SpatialPGM<int> index(/*epsilon_x=*/1, /*epsilon_y=*/1);
  index.build(points);

  // Check that x-index was built (has segments)
  EXPECT_GT(index.x_index().segments_count(), 0u);
  
  // Check that y-indexes were built for x segments
  const auto& x_segments = index.x_segments();
  ASSERT_GT(x_segments.size(), 0u);
  
  // Each x segment should have a y-index
  for (const auto& x_seg : x_segments) {
    EXPECT_GT(x_seg.y_index.segments_count(), 0u);
    EXPECT_FALSE(x_seg.points.empty());
  }

  // Check that we have the expected number of segments
  EXPECT_GT(index.num_x_segments(), 0u);
  
  // Check that points are properly distributed
  if (x_segments.size() > 0) {
    EXPECT_GT(index.num_y_segments(0), 0u);
  }

  // Verify the index contains the points
  for (const auto& point : points) {
    EXPECT_TRUE(index.exists(point));
  }
}

TEST(SpatialPGMTest, RangeQueryFiltersByAxis) {
  std::vector<Point<int>> points = {
      {0, 0}, {1, 0}, {2, 0}, {10, 1}, {11, 2}, {12, 2}};

  SpatialPGM<int> index(/*epsilon_x=*/1, /*epsilon_y=*/1);
  index.build(points);

  auto result_one = index.range_query({0, 0}, {11, 1});
  sort_points(result_one);
  std::vector<Point<int>> expected_one = {{0, 0}, {1, 0}, {2, 0}, {10, 1}};
  ASSERT_EQ(result_one.size(), expected_one.size());
  for (size_t i = 0; i < result_one.size(); ++i) {
    EXPECT_EQ(result_one[i].x, expected_one[i].x);
    EXPECT_EQ(result_one[i].y, expected_one[i].y);
  }

  auto result_two = index.range_query({10, 2}, {12, 2});
  sort_points(result_two);
  std::vector<Point<int>> expected_two = {{11, 2}, {12, 2}};
  ASSERT_EQ(result_two.size(), expected_two.size());
  for (size_t i = 0; i < result_two.size(); ++i) {
    EXPECT_EQ(result_two[i].x, expected_two[i].x);
    EXPECT_EQ(result_two[i].y, expected_two[i].y);
  }

  auto result_none = index.range_query({20, 0}, {30, 5});
  EXPECT_TRUE(result_none.empty());
}

TEST(SpatialPGMTest, ExistsFunction) {
  std::vector<Point<int>> points = {
      {0, 0}, {1, 0}, {2, 0}, {10, 1}, {11, 2}, {12, 2}};

  SpatialPGM<int> index(/*epsilon_x=*/1, /*epsilon_y=*/1);
  index.build(points);

  // Test that all built points exist
  for (const auto& point : points) {
    EXPECT_TRUE(index.exists(point)) 
        << "Point (" << point.x << ", " << point.y << ") should exist";
  }

  // Test that non-existent points don't exist
  EXPECT_FALSE(index.exists({0, 1}));
  EXPECT_FALSE(index.exists({1, 1}));
  EXPECT_FALSE(index.exists({5, 0}));
  EXPECT_FALSE(index.exists({10, 0}));
  EXPECT_FALSE(index.exists({11, 1}));
  EXPECT_FALSE(index.exists({12, 1}));
  EXPECT_FALSE(index.exists({20, 0}));
  EXPECT_FALSE(index.exists({0, 10}));
}

TEST(SpatialPGMTest, ExistsWithEmptyIndex) {
  SpatialPGM<int> index(/*epsilon_x=*/1, /*epsilon_y=*/1);
  std::vector<Point<int>> empty_points;
  index.build(empty_points);

  // Test that exists returns false for empty index
  EXPECT_FALSE(index.exists({0, 0}));
  EXPECT_FALSE(index.exists({1, 1}));
}

TEST(SpatialPGMTest, ExistsWithSinglePoint) {
  std::vector<Point<int>> points = {{5, 5}};

  SpatialPGM<int> index(/*epsilon_x=*/1, /*epsilon_y=*/1);
  index.build(points);

  // Test that the single point exists
  EXPECT_TRUE(index.exists({5, 5}));

  // Test that other points don't exist
  EXPECT_FALSE(index.exists({5, 4}));
  EXPECT_FALSE(index.exists({5, 6}));
  EXPECT_FALSE(index.exists({4, 5}));
  EXPECT_FALSE(index.exists({6, 5}));
}

TEST(SpatialPGMTest, ExistsWithDuplicatePoints) {
  std::vector<Point<int>> points = {
      {0, 0}, {0, 0}, {1, 1}, {1, 1}, {2, 2}};

  SpatialPGM<int> index(/*epsilon_x=*/1, /*epsilon_y=*/1);
  index.build(points);

  // Test that duplicate points exist
  EXPECT_TRUE(index.exists({0, 0}));
  EXPECT_TRUE(index.exists({1, 1}));
  EXPECT_TRUE(index.exists({2, 2}));

  // Test that non-existent points don't exist
  EXPECT_FALSE(index.exists({0, 1}));
  EXPECT_FALSE(index.exists({3, 3}));
}

TEST(SpatialPGMTest, ExistsLargeIndex) {
  int num_points = 100000;
  std::vector<Point<int>> points;
  for (int i = 0; i < num_points; ++i) {
    points.push_back({rand(), rand()});
  }
  SpatialPGM<int> index(/*epsilon_x=*/1, /*epsilon_y=*/1);
  index.build(points);
  for (const auto& point : points) {
    EXPECT_TRUE(index.exists(point));
  }
}

// Helper function to manually find points in range
template <class T>
std::vector<Point<T>> manual_range_query(const std::vector<Point<T>>& points,
                                          const Point<T>& lower,
                                          const Point<T>& upper) {
  std::vector<Point<T>> result;
  for (const auto& p : points) {
    if (p.x >= lower.x && p.x <= upper.x && p.y >= lower.y && p.y <= upper.y) {
      result.push_back(p);
    }
  }
  sort_points(result);
  return result;
}

TEST(SpatialPGMTest, RangeQueryWithRandomPoints) {
  const int num_points = 10000;
  const int num_queries = 100;
  
  // Generate random points
  std::vector<Point<int>> points;
  int min_x = INT_MAX, min_y = INT_MAX, max_x = INT_MIN, max_y = INT_MIN;
  
  srand(42);  // Fixed seed for reproducibility
  for (int i = 0; i < num_points; ++i) {
    points.push_back((Point<int>){rand(), rand()});
    min_x = std::min(min_x, points[i].x);
    min_y = std::min(min_y, points[i].y);
    max_x = std::max(max_x, points[i].x);
    max_y = std::max(max_y, points[i].y);
  }
  
  // Remove duplicates
  sort_points(points);
  points.erase(std::unique(points.begin(), points.end()), points.end());
  
  // Build index
  SpatialPGM<int> index(/*epsilon_x=*/64, /*epsilon_y=*/64);
  index.build(points);
  
  // Test full range query
  auto full_result = index.range_query({min_x, min_y}, {max_x, max_y});
  auto expected_full = manual_range_query(points, {min_x, min_y}, {max_x, max_y});
  sort_points(full_result);
  
  ASSERT_EQ(full_result.size(), expected_full.size())
      << "Full range query should return all points";
  for (size_t i = 0; i < full_result.size(); ++i) {
    EXPECT_EQ(full_result[i].x, expected_full[i].x);
    EXPECT_EQ(full_result[i].y, expected_full[i].y);
  }
  
  // Test multiple random query rectangles
  srand(123);  // Fixed seed for reproducibility
  // Generate rectangle query
  std::vector<std::pair<Point<int>, Point<int>>> query_data;
  for (int i = 0; i < num_queries; ++i) {
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
  for (int q = 0; q < num_queries; ++q) {
    auto [lower, upper] = query_data[q];

    // Query using index
    auto index_result = index.range_query(lower, upper);
    sort_points(index_result);
    
    // Query using manual inspection
    auto manual_result = manual_range_query(points, lower, upper);
    
    // Compare results
    ASSERT_EQ(index_result.size(), manual_result.size());
    
    for (size_t i = 0; i < index_result.size(); ++i) {
      EXPECT_EQ(index_result[i].x, manual_result[i].x);
      EXPECT_EQ(index_result[i].y, manual_result[i].y);
    }
  }
}

TEST(SpatialPGMTest, RangeQueryWithEmptyResults) {
  const int num_points = 1000;
  
  // Generate random points in a specific range
  std::vector<Point<int>> points;
  srand(456);  // Fixed seed for reproducibility
  for (int i = 0; i < num_points; ++i) {
    int x = 1000 + (rand() % 1000);  // x in [1000, 2000)
    int y = 1000 + (rand() % 1000);  // y in [1000, 2000)
    points.push_back({x, y});
  }
  
  // Build index
  SpatialPGM<int> index(/*epsilon_x=*/64, /*epsilon_y=*/64);
  index.build(points);
  
  // Test queries that should return empty results
  auto result1 = index.range_query({0, 0}, {500, 500});
  EXPECT_TRUE(result1.empty()) << "Query outside point range should return empty";
  
  auto result2 = index.range_query({3000, 3000}, {4000, 4000});
  EXPECT_TRUE(result2.empty()) << "Query outside point range should return empty";
  
  auto result3 = index.range_query({500, 500}, {999, 999});
  EXPECT_TRUE(result3.empty()) << "Query just before point range should return empty";
  
  auto result4 = index.range_query({2001, 2001}, {3000, 3000});
  EXPECT_TRUE(result4.empty()) << "Query just after point range should return empty";
}
