#ifndef SPATIALPGM_H
#define SPATIALPGM_H

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "point.h"
#include "pgm/pgm_index.hpp"

// A spatial variant of the PGM-Index. It first builds a hierarchical PGM index
// on the x-axis (PGMx). For every leaf x-segment it then builds another
// PGM index on the y values of the points belonging to that x segment (PGMaY).
template <class T, size_t EpsilonX = 64, size_t EpsilonRecursiveX = 4,
          size_t EpsilonY = 64, size_t EpsilonRecursiveY = 4>
class SpatialPGM {
 public:
  // Represents a leaf segment in the x-index along with its y-index
  struct XSegment {
    size_t start;  // Start index in points_by_x_
    size_t end;    // End index in points_by_x_
    pgm::PGMIndex<T, EpsilonY, EpsilonRecursiveY> y_index;  // PGM index over y
    std::vector<Point<T>> points;  // Points sorted by y for this segment
  };

  SpatialPGM() = default;

  explicit SpatialPGM(size_t epsilon_x, size_t epsilon_y)
      : epsilon_x_(epsilon_x), epsilon_y_(epsilon_y) {}

  void build(const std::vector<Point<T>>& points) {
    points_by_x_ = points;
    x_segments_.clear();

    if (points_by_x_.empty()) {
      return;
    }

    // Sort points by x (and by y for ties)
    std::sort(points_by_x_.begin(), points_by_x_.end(),
              [](const Point<T>& lhs, const Point<T>& rhs) {
                if (lhs.x == rhs.x) return lhs.y < rhs.y;
                return lhs.x < rhs.x;
              });

    // Extract x coordinates for building the x-index
    std::vector<T> xs;
    xs.reserve(points_by_x_.size());
    for (const auto& p : points_by_x_) {
      xs.push_back(p.x);
    }

    // Build hierarchical PGM index over x coordinates
    if constexpr (EpsilonX != 64 || EpsilonRecursiveX != 4) {
      x_index_ = pgm::PGMIndex<T, EpsilonX, EpsilonRecursiveX>(xs.begin(), xs.end());
    } else {
      x_index_ = pgm::PGMIndex<T, 64, 4>(xs.begin(), xs.end());
    }

    // Get number of leaf segments
    size_t segments_count = x_index_.segments_count();
    if (segments_count == 0) {
      return;
    }

    // Map points to segments based on the PGM index's actual segment boundaries.
    // The PGM index segments have different sizes based on the key distribution,
    // so we cannot simply divide points evenly. We use the actual segment key
    // boundaries to determine which points belong to which segment.
    
    // Get the key boundaries for all leaf segments
    std::vector<T> segment_keys = x_index_.get_leaf_segment_keys();
    if (segment_keys.size() != segments_count) {
      // Fallback if something went wrong
      segment_keys.clear();
      fprintf(stderr, "Number of segments in x is incorrect\n");
      for (size_t i = 0; i < segments_count; ++i) {
        segment_keys.push_back(x_index_.get_leaf_segment_key(i));
      }
    }
    
    std::vector<std::vector<Point<T>>> segment_points(segments_count);
    // Track the start and end indices in points_by_x_ for each segment
    std::vector<size_t> segment_starts(segments_count);
    std::vector<size_t> segment_ends(segments_count);
    
    // Since both points_by_x_ and segment_keys are sorted, we can process them
    // in a single pass using a two-pointer approach. This is O(n + m) where
    // n is the number of points and m is the number of segments.
    //
    // The algorithm: iterate through segments (outer loop), and for each segment,
    // iterate through points (inner loop) assigning them to the current segment
    // until we reach a point that belongs to the next segment.
    size_t point_idx = 0;
    for (size_t seg_idx = 0; seg_idx < segments_count; ++seg_idx) {
      // Determine the key boundary for the next segment (or sentinel for last segment)
      T next_segment_key = (seg_idx + 1 < segment_keys.size()) 
                          ? segment_keys[seg_idx + 1] 
                          : std::numeric_limits<T>::max();
      
      // Record the start index for this segment
      segment_starts[seg_idx] = point_idx;
      
      // Assign points to this segment until we reach a point that belongs to the next segment
      while (point_idx < points_by_x_.size() && 
             points_by_x_[point_idx].x < next_segment_key) {
        segment_points[seg_idx].push_back(points_by_x_[point_idx]);
        ++point_idx;
      }
      
      // Record the end index (one past the last point assigned to this segment)
      segment_ends[seg_idx] = point_idx;
    }

    x_segments_.reserve(segments_count);
    
    // Build y-indexes for each segment
    for (size_t seg_idx = 0; seg_idx < segments_count; ++seg_idx) {
      if (segment_points[seg_idx].empty()) {
        continue;
      }

      // Extract points for this segment
      std::vector<Point<T>> segment_pts = std::move(segment_points[seg_idx]);

      if (segment_pts.empty()) {
        continue;
      }

      // Sort points by y for this segment
      std::sort(segment_pts.begin(), segment_pts.end(),
                [](const Point<T>& lhs, const Point<T>& rhs) {
                  if (lhs.y == rhs.y) return lhs.x < rhs.x;
                  return lhs.y < rhs.y;
                });

      // Extract y coordinates
      std::vector<T> ys;
      ys.reserve(segment_pts.size());
      for (const auto& p : segment_pts) {
        ys.push_back(p.y);
      }

      if (ys.empty()) {
        continue;
      }

      // Build PGM index over y coordinates
      XSegment x_seg;
      // Use the start and end indices we tracked during the mapping phase
      x_seg.start = segment_starts[seg_idx];
      x_seg.end = segment_ends[seg_idx];
      x_seg.points = std::move(segment_pts);
      
      if constexpr (EpsilonY != 64 || EpsilonRecursiveY != 4) {
        x_seg.y_index = pgm::PGMIndex<T, EpsilonY, EpsilonRecursiveY>(ys.begin(), ys.end());
      } else {
        x_seg.y_index = pgm::PGMIndex<T, 64, 4>(ys.begin(), ys.end());
      }

      x_segments_.push_back(std::move(x_seg));
    }
  }

  size_t num_x_segments() const {
    return x_segments_.size();
  }

  size_t num_y_segments(size_t x_segment_index) const {
    if (x_segment_index >= x_segments_.size()) {
      throw std::out_of_range("x_segment_index out of range");
    }
    return x_segments_[x_segment_index].y_index.segments_count();
  }

  const std::vector<XSegment>& x_segments() const { return x_segments_; }
  const pgm::PGMIndex<T, EpsilonX, EpsilonRecursiveX>& x_index() const { return x_index_; }

  // Check if a point exists in the index
  // Similar to PGMIndex::search, first lookup in x_index_, then in the corresponding y_index
  bool exists(const Point<T>& key) const {
    if (x_segments_.empty() || points_by_x_.empty()) {
      return false;
    }

    // Step 1: Use x_index_ to find which segment contains this point's x coordinate
    auto x_approx = x_index_.search(key.x);
    
    // Search in the range [x_approx.lo, x_approx.hi) to find the segment
    size_t search_start = std::min(x_approx.lo, points_by_x_.size());
    size_t search_end = std::min(x_approx.hi, points_by_x_.size());
    
    // Find the segment that contains this point's x coordinate
    for (const auto& x_seg : x_segments_) {
      // Check if this segment overlaps with the search range
      if (x_seg.end <= search_start || x_seg.start >= search_end) {
        continue;
      }

      // Verify the point's x is actually in this segment's range
      // by checking the actual x values in the segment
      bool x_found = false;
      for (size_t i = x_seg.start; i < x_seg.end && i < points_by_x_.size(); ++i) {
        if (points_by_x_[i].x == key.x) {
          x_found = true;
          break;
        }
        if (points_by_x_[i].x > key.x) {
          break;  // Points are sorted by x, so we can stop
        }
      }
      
      if (!x_found) {
        continue;
      }

      // Step 2: Use the y_index for this segment to find the point
      auto y_approx = x_seg.y_index.search(key.y);
      
      // Search in the y range [y_approx.lo, y_approx.hi)
      size_t y_search_start = std::min(y_approx.lo, x_seg.points.size());
      size_t y_search_end = std::min(y_approx.hi, x_seg.points.size());
      
      // Perform binary search to find the exact point, similar to PGMIndex
      auto begin_it = x_seg.points.begin() + static_cast<std::ptrdiff_t>(y_search_start);
      auto end_it = x_seg.points.begin() + static_cast<std::ptrdiff_t>(y_search_end);
      
      // Use lower_bound to find the first point with y >= key.y
      auto point_it = std::lower_bound(begin_it, end_it, key.y,
                                       [](const Point<T>& p, T y) { return p.y < y; });
      
      // Check if the exact point exists
      for (; point_it != end_it && point_it->y == key.y; ++point_it) {
        if (point_it->x == key.x) {
          return true;
        }
      }
    }

    return false;
  }

  std::vector<Point<T>> range_query(const Point<T>& lower,
                                    const Point<T>& upper) const {
    std::vector<Point<T>> results;
    
    if (x_segments_.empty() || points_by_x_.empty()) {
      return results;
    }

    // For each x segment, check if it overlaps with the x range
    for (const auto& x_seg : x_segments_) {
      // Optimization 1: Skip segments that are completely outside the query range
      // Since points_by_x_ is sorted by x, we can check segment boundaries directly
      if (x_seg.start >= points_by_x_.size() || x_seg.end == 0) {
        continue;
      }
      
      // Check if segment is completely before or after the query range
      // This avoids iterating through all points in the segment
      T seg_min_x = points_by_x_[x_seg.start].x;
      T seg_max_x = points_by_x_[x_seg.end - 1].x;
      
      if (seg_max_x < lower.x || seg_min_x > upper.x) {
        continue;  // Segment is completely outside the query range, skip it
      }

      // Optimization 2: Check if the entire segment is within the x range
      // If so, we don't need to check x when iterating through y-sorted points
      bool segment_fully_in_x_range = (seg_min_x >= lower.x && seg_max_x <= upper.x);

      // Use y-index to find points in y range
      auto y_approx = x_seg.y_index.search(lower.y);
      
      // Search in the y range [y_approx.lo, y_approx.hi)
      size_t y_search_start = std::min(y_approx.lo, x_seg.points.size());
      size_t y_search_end = std::min(y_approx.hi, x_seg.points.size());
      
      // Perform binary search to find the actual y range
      auto y_begin_it = x_seg.points.begin() + static_cast<std::ptrdiff_t>(y_search_start);
      auto y_end_it = x_seg.points.begin() + static_cast<std::ptrdiff_t>(y_search_end);
      
      auto it = std::lower_bound(y_begin_it, y_end_it, lower.y,
                                 [](const Point<T>& p, T y) { return p.y < y; });
      
      // Collect all points in the range
      if (segment_fully_in_x_range) {
        // No need to check x - all points in this segment are in the x range
        for (; it != x_seg.points.end() && it->y <= upper.y; ++it) {
          results.push_back(*it);
        }
      } else {
        // Need to check x because segment partially overlaps with x range
        for (; it != x_seg.points.end() && it->y <= upper.y; ++it) {
          if (it->x >= lower.x && it->x <= upper.x) {
            results.push_back(*it);
          }
        }
      }
    }

    return results;
  }

 private:

  size_t epsilon_x_ = EpsilonX;
  size_t epsilon_y_ = EpsilonY;
  std::vector<Point<T>> points_by_x_;
  pgm::PGMIndex<T, EpsilonX, EpsilonRecursiveX> x_index_;
  std::vector<XSegment> x_segments_;
};

#endif  // SPATIALPGM_H
