#ifndef __MPGM_H__
#define __MPGM_H__

#include <algorithm>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <stack>
#include <type_traits>
#include <utility>
#include <vector>

#include "point.h"
#include "pgm/piecewise_linear_model.hpp"

// Alternating-axis spatial index built from stacked PGM indexes.
// Level 0 uses x as the key; each leaf segment spawns a child built on y;
// the next level flips back to x, and so on.
template <class T, size_t EpsilonX = 64, size_t EpsilonRecursiveX = 4,
          size_t EpsilonY = 64, size_t EpsilonRecursiveY = 4>
class mPGM {
  enum class Axis { X = 0, Y = 1 };

  // Segment struct similar to PGM's Segment, storing key, slope, and intercept
  struct Segment {
    T key;              ///< The first key that the segment indexes.
    double slope;        ///< The slope of the segment.
    double intercept;    ///< The intercept of the segment.

    Segment() = default;
    Segment(T key, double slope, double intercept) : key(key), slope(slope), intercept(intercept) {}

    explicit Segment(const typename pgm::internal::OptimalPiecewiseLinearModel<T, size_t>::CanonicalSegment &cs)
        : key(cs.get_first_x()) {
      auto [cs_slope, cs_intercept] = cs.get_floating_point_segment(key);
      slope = static_cast<double>(cs_slope);
      intercept = static_cast<double>(cs_intercept);
    }

    friend inline bool operator<(const Segment &s, const T &k) { return s.key < k; }
    friend inline bool operator<(const T &k, const Segment &s) { return k < s.key; }
    friend inline bool operator<(const Segment &s, const Segment &t) { return s.key < t.key; }

    operator T() { return key; }

    /**
     * Returns the approximate position of the specified key.
     * This follows the same approach as PGM's Segment::operator().
     */
    inline size_t operator()(const T &k) const {
      size_t pos;
      if constexpr (std::is_same_v<T, int64_t> || std::is_same_v<T, int32_t>)
        pos = size_t(slope * double(std::make_unsigned_t<T>(k) - key));
      else
        pos = size_t(slope * double(k - key));
      return pos + static_cast<size_t>(intercept);
    }
  };

  struct Node {
    Axis axis;
    bool is_leaf = false;
    std::vector<Point<T>> points;  // Points sorted by current axis

    // Segment metadata
    std::vector<Segment> segments;  // PGM segments with key, slope, and intercept
    std::vector<size_t> segment_start;  // indices into points (sorted by axis)
    std::vector<size_t> segment_end;    // indices into points (sorted by axis)
    std::vector<size_t> segment_upper;  // exclusive upper bound for segment positions
    std::vector<T> segment_min;       // min key in segment along current axis
    std::vector<T> segment_max;       // max key in segment along current axis
    std::vector<Node*> children;      // child for each segment (nullptr for leaf)
  };

 public:
  mPGM() = default;
  ~mPGM() { destroy(root_); }

  void bulk_load(Point<T>* vs, size_t num_keys) {
    std::vector<Point<T>> pts(vs, vs + num_keys);
    build(pts);
  }

  void build(const std::vector<Point<T>>& points) {
    destroy(root_);
    data_ = points;
    root_ = build_node(data_, Axis::X);
  }

  void insert(const Point<T>& point) {
    data_.push_back(point);
    destroy(root_);
    root_ = build_node(data_, Axis::X);
  }

  size_t index_size() const {
    size_t bytes = sizeof(data_) + data_.capacity() * sizeof(Point<T>);
    if (!root_) return bytes;

    std::stack<const Node*> s;
    s.push(root_);
    while (!s.empty()) {
      const Node* node = s.top();
      s.pop();

      bytes += sizeof(*node);
      bytes += node->points.capacity() * sizeof(Point<T>);
      bytes += node->segments.capacity() * sizeof(Segment);
      bytes += node->segment_start.capacity() * sizeof(size_t);
      bytes += node->segment_end.capacity() * sizeof(size_t);
      bytes += node->segment_upper.capacity() * sizeof(size_t);
      bytes += node->segment_min.capacity() * sizeof(T);
      bytes += node->segment_max.capacity() * sizeof(T);
      bytes += node->children.capacity() * sizeof(Node*);

      for (const auto* child : node->children) {
        if (child) s.push(child);
      }
    }
    return bytes;
  }

  size_t size() const { return data_.size(); }
  bool empty() const { return data_.empty(); }

  bool exists(const Point<T>& key) const {
    if (!root_) return false;

    // If root is a leaf, we cannot rely on a model; defer to leaf scan.
    if (root_->is_leaf) {
      return exists_node(root_, key, kNoPrediction);
    }

    // Seed the lookup with the root's linear model prediction so every level
    // is guided by its corresponding segment.
    T axis_key = axis_value(key, root_->axis);
    size_t seg_idx = segment_for_key(root_->segments, axis_key);
    size_t predicted = root_->segments[seg_idx](axis_key);
    size_t upper_exclusive =
        (seg_idx + 1 < root_->segment_upper.size()) ? root_->segment_upper[seg_idx]
                                                    : root_->points.size();
    if (upper_exclusive > 0) {
      predicted = std::min(predicted, upper_exclusive - 1);
    }
    if (!root_->points.empty()) {
      predicted = std::min(predicted, root_->points.size() - 1);
    }
    return exists_node(root_, key, predicted);
  }

  void rangeQuery(const Point<T>& lower, const Point<T>& upper,
                  std::vector<Point<T>>& out) const {
    range_query_node(root_, lower, upper, out);
  }

  std::vector<Point<T>> range_query(const Point<T>& lower,
                                    const Point<T>& upper) const {
    std::vector<Point<T>> res;
    range_query_node(root_, lower, upper, res);
    return res;
  }

  int range_query(const Point<T>& lower, const Point<T>& upper,
                  Point<T>* results) const {
    auto vec = range_query(lower, upper);
    std::copy(vec.begin(), vec.end(), results);
    return static_cast<int>(vec.size());
  }

 private:
  static constexpr size_t LEAF_SIZE = 32;

  static T axis_value(const Point<T>& p, Axis axis) {
    return axis == Axis::X ? p.x : p.y;
  }

  Node* build_node(std::vector<Point<T>> pts, Axis axis) {
    if (pts.empty()) {
      return nullptr;
    }

    // Sort points by current axis, then by the other axis for stability.
    std::sort(pts.begin(), pts.end(),
              [axis](const Point<T>& a, const Point<T>& b) {
                if (axis == Axis::X) {
                  if (a.x == b.x) return a.y < b.y;
                  return a.x < b.x;
                } else {
                  if (a.y == b.y) return a.x < b.x;
                  return a.y < b.y;
                }
              });

    Node* node = new Node();
    node->axis = axis;
    node->points = std::move(pts);

    // Use PGM's segmentation logic to create segments along current axis
    // Even for leaf nodes, we create segments to use linear models for prediction
    size_t n = node->points.size();
    size_t epsilon = (axis == Axis::X) ? EpsilonX : EpsilonY;
    
    // Mark as leaf if size is small enough
    if (node->points.size() <= LEAF_SIZE) {
      node->is_leaf = true;
    }
    
    // Input function: get key at index i
    auto in_fun = [&node, axis](size_t i) -> T {
      return axis_value(node->points[i], axis);
    };
    
    // Output function: create Segment objects directly
    std::vector<Segment> temp_segments;
    auto out_fun = [&temp_segments](const auto& cs) {
      temp_segments.emplace_back(cs);
    };
    
    // Create segments using PGM's segmentation logic
    size_t segments_count = pgm::internal::make_segmentation(n, epsilon, in_fun, out_fun);
    
    if (segments_count == 0) {
      node->is_leaf = true;
      return node;
    }
    
    node->segments = std::move(temp_segments);

    node->segment_start.resize(segments_count);
    node->segment_end.resize(segments_count);
    node->segment_min.resize(segments_count);
    node->segment_max.resize(segments_count);
    node->children.resize(segments_count, nullptr);
    node->segment_upper.resize(segments_count);

    std::vector<std::vector<Point<T>>> segment_points(segments_count);
    size_t point_idx = 0;
    for (size_t seg_idx = 0; seg_idx < segments_count; ++seg_idx) {
      T next_key = (seg_idx + 1 < node->segments.size())
                       ? node->segments[seg_idx + 1].key
                       : std::numeric_limits<T>::max();
      node->segment_start[seg_idx] = point_idx;
      while (point_idx < node->points.size() &&
             axis_value(node->points[point_idx], axis) < next_key) {
        segment_points[seg_idx].push_back(node->points[point_idx]);
        ++point_idx;
      }
      node->segment_end[seg_idx] = point_idx;
      if (!segment_points[seg_idx].empty()) {
        node->segment_min[seg_idx] =
            axis_value(segment_points[seg_idx].front(), axis);
        node->segment_max[seg_idx] =
            axis_value(segment_points[seg_idx].back(), axis);
      } else {
        node->segment_min[seg_idx] = std::numeric_limits<T>::max();
        node->segment_max[seg_idx] = std::numeric_limits<T>::lowest();
      }
      // Upper bound (exclusive) for positions in this segment matches the start
      // of the next segment, or total size for the last segment.
      node->segment_upper[seg_idx] =
          (seg_idx + 1 < segments_count) ? point_idx : node->points.size();
    }

    size_t non_empty_segments = 0;
    for (const auto& segment : segment_points) {
      if (!segment.empty()) ++non_empty_segments;
    }

    // If segmentation failed to split the data, fall back to a simple median
    // split to guarantee progress and avoid infinite recursion.
    if (non_empty_segments <= 1) {
      if (node->is_leaf || node->points.size() <= LEAF_SIZE) {
        node->is_leaf = true;
        return node;
      }

      const size_t mid = node->points.size() / 2;
      std::vector<Point<T>> left(node->points.begin(),
                                 node->points.begin() + mid);
      std::vector<Point<T>> right(node->points.begin() + mid,
                                  node->points.end());

      // For fallback median split, compute simple linear approximations
      // Slope = (end_index - start_index) / (max_key - min_key)
      T left_min = axis_value(left.front(), axis);
      T left_max = axis_value(left.back(), axis);
      T right_min = axis_value(right.front(), axis);
      T right_max = axis_value(right.back(), axis);
      double left_slope = (left_max != left_min) ? static_cast<double>(left.size() - 1) / static_cast<double>(left_max - left_min) : 0.0;
      double right_slope = (right_max != right_min) ? static_cast<double>(right.size() - 1) / static_cast<double>(right_max - right_min) : 0.0;
      double left_intercept = 0.0;  // intercept at left_min
      double right_intercept = static_cast<double>(mid);  // intercept at right_min
      node->segments = {
        Segment(left_min, left_slope, left_intercept),
        Segment(right_min, right_slope, right_intercept)
      };
      node->segment_start = {0, mid};
      node->segment_end = {mid, node->points.size()};
      node->segment_upper = {mid, node->points.size()};
      node->segment_min = {left_min, right_min};
      node->segment_max = {left_max, right_max};
      node->children.assign(2, nullptr);

      Axis child_axis = (axis == Axis::X) ? Axis::Y : Axis::X;
      node->children[0] = build_node(std::move(left), child_axis);
      node->children[1] = build_node(std::move(right), child_axis);
      return node;
    }

    // Only create children if this is not a leaf node
    if (!node->is_leaf) {
      Axis child_axis = (axis == Axis::X) ? Axis::Y : Axis::X;
      for (size_t seg_idx = 0; seg_idx < segments_count; ++seg_idx) {
        if (!segment_points[seg_idx].empty()) {
          node->children[seg_idx] =
              build_node(std::move(segment_points[seg_idx]), child_axis);
        }
      }
    }

    return node;
  }

  static size_t segment_for_key(const std::vector<Segment>& segments, T key) {
    // Find the rightmost segment with key <= the sought key
    auto it = std::upper_bound(segments.begin(), segments.end(), key);
    if (it == segments.begin()) {
      return 0;
    }
    return static_cast<size_t>(std::distance(segments.begin(), it - 1));
  }

  static size_t segment_for_position(const Node* node, size_t pos) {
    auto it =
        std::upper_bound(node->segment_start.begin(), node->segment_start.end(), pos);
    if (it == node->segment_start.begin()) {
      return 0;
    }
    return static_cast<size_t>(std::distance(node->segment_start.begin(), it - 1));
  }

  static constexpr size_t kNoPrediction = std::numeric_limits<size_t>::max();

  bool exists_node(const Node* node, const Point<T>& key,
                   size_t predicted_pos = kNoPrediction) const {
    if (!node) return false;

    if (node->is_leaf) {
      if (node->points.empty()) return false;

      // Use linear model from segments to predict position, just like PGM does
      size_t epsilon = (node->axis == Axis::X) ? EpsilonX : EpsilonY;
      T axis_key = axis_value(key, node->axis);
      
      size_t center = predicted_pos;
      if (center == kNoPrediction) {
        // If no prediction from parent, use segment's linear model
        if (!node->segments.empty()) {
          size_t seg_idx = segment_for_key(node->segments, axis_key);
          center = node->segments[seg_idx](axis_key);
          // Clamp using segment boundary (similar to PGM)
          if (seg_idx + 1 < node->segments.size()) {
            center = std::min(center, static_cast<size_t>(node->segments[seg_idx + 1].intercept));
          }
          center = std::min(center, node->points.size() - 1);
        } else {
          // Fallback to binary search if no segments (shouldn't happen normally)
          auto it = std::lower_bound(
              node->points.begin(), node->points.end(), axis_key,
              [axis = node->axis](const Point<T>& p, const T& v) {
                return axis_value(p, axis) < v;
              });
          center = static_cast<size_t>(std::distance(node->points.begin(), it));
          if (center == node->points.size()) {
            center = node->points.size() - 1;
          }
        }
      } else {
        center = std::min(center, node->points.size() - 1);
      }

      // Create search window using epsilon, matching PGM's approach
      size_t lo = (center > epsilon) ? center - epsilon : 0;
      size_t hi = std::min(center + epsilon + 2, node->points.size());
      for (size_t i = lo; i < hi; ++i) {
        const auto& p = node->points[i];
        if (p.x == key.x && p.y == key.y) {
          return true;
        }
      }
      return false;
    }

    // Predict the target position with the current level's linear model.
    T axis_key = axis_value(key, node->axis);
    size_t responsible_seg = segment_for_key(node->segments, axis_key);
    size_t predicted_pos_local = node->segments[responsible_seg](axis_key);
    // Clamp prediction using the next segment's boundary, mirroring PGM search.
    size_t upper_exclusive = node->segment_upper[responsible_seg];
    if (upper_exclusive == 0) return false;
    predicted_pos_local = std::min(predicted_pos_local, upper_exclusive - 1);
    if (!node->points.empty()) {
      predicted_pos_local = std::min(predicted_pos_local, node->points.size() - 1);
    }

    // Map predicted position to segment index.
    size_t predicted_seg = segment_for_position(node, predicted_pos_local);

    auto descend = [&](size_t seg_idx) -> bool {
      if (seg_idx >= node->children.size() || !node->children[seg_idx]) return false;
      if (node->segment_min[seg_idx] > axis_key ||
          axis_key > node->segment_max[seg_idx]) {
        return false;
      }
      size_t child_predicted =
          (predicted_pos_local > node->segment_start[seg_idx])
              ? predicted_pos_local - node->segment_start[seg_idx]
              : 0;
      return exists_node(node->children[seg_idx], key, child_predicted);
    };

    // Try the predicted segment first, then immediate neighbors if needed.
    if (descend(predicted_seg)) return true;
    if (predicted_seg > 0 && descend(predicted_seg - 1)) return true;
    if (predicted_seg + 1 < node->segments.size() && descend(predicted_seg + 1))
      return true;

    return false;
  }

  void range_query_node(const Node* node, const Point<T>& lower,
                        const Point<T>& upper,
                        std::vector<Point<T>>& out) const {
    if (!node) return;
    if (node->is_leaf) {
      for (const auto& p : node->points) {
        if (p.x >= lower.x && p.x <= upper.x &&
            p.y >= lower.y && p.y <= upper.y) {
          out.push_back(p);
        }
      }
      return;
    }

    Axis axis = node->axis;
    T q_lower = (axis == Axis::X) ? lower.x : lower.y;
    T q_upper = (axis == Axis::X) ? upper.x : upper.y;

    // Find the range of segments that might overlap with the query range
    // Find the first segment that might contain q_lower
    auto lower_it = std::lower_bound(node->segments.begin(), 
                                     node->segments.end(), q_lower);
    size_t start_seg = 0;
    if (lower_it != node->segments.begin()) {
      // Include the segment before if it might overlap
      start_seg = static_cast<size_t>(
          std::distance(node->segments.begin(), lower_it) - 1);
    }
    
    // Find the last segment that might contain q_upper
    auto upper_it = std::upper_bound(node->segments.begin(), 
                                      node->segments.end(), q_upper);
    size_t end_seg = static_cast<size_t>(
        std::distance(node->segments.begin(), upper_it));

    // Check all segments in the computed range
    for (size_t i = start_seg; i < end_seg && i < node->children.size(); ++i) {
      if (!node->children[i]) continue;
      // Check if segment overlaps with query range on current axis
      if (node->segment_max[i] < q_lower || node->segment_min[i] > q_upper) {
        continue;
      }
      range_query_node(node->children[i], lower, upper, out);
    }
  }

  void destroy(Node* node) {
    if (!node) return;
    for (auto* child : node->children) {
      destroy(child);
    }
    delete node;
  }

  Node* root_ = nullptr;
  std::vector<Point<T>> data_;
};

#endif // __MPGM_H__
