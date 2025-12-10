#ifndef BTREE_ZORDER_HPP
#define BTREE_ZORDER_HPP

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

#include "btree_olc.h"
#include "morton.h"
#include "point.h"

using namespace libmorton;
using namespace btree_olc_internal;

// A simple B+-tree based index over 2D points using Z-order (Morton) encoding.
// Points are encoded into a 64-bit Morton key and indexed with the BTree
// implementation from btree_olc.h. 
//
// Note: Multiple points can encode to the same z-order number:
//   - For integer coordinates: Points with identical (x, y) values
//   - For double coordinates: Points that map to the same uint32_t after
//     precision multiplication (precision = 2^32 - 1)
// For duplicate Morton keys, we keep a list of the original points in payload
// buckets and filter during lookups/range queries to find exact matches.
template <class T>
class BTreeZOrder {
 public:
 BTreeZOrder() : btree_() {}

  ~BTreeZOrder() {
    destroy_tree(btree_.root);
    btree_.root = nullptr;
  }

  void set_bounds(double min_x, double min_y, double max_x, double max_y) {
    bounds_.min_x = min_x;
    bounds_.min_y = min_y;
    bounds_.span_x = std::max(max_x - min_x, std::numeric_limits<double>::epsilon());
    bounds_.span_y = std::max(max_y - min_y, std::numeric_limits<double>::epsilon());
  }

  void bulk_load(Point<T>* vs, int num_keys) {
    reset();
    if (num_keys <= 0) {
      return;
    }

    std::vector<std::pair<ZKey, Point<T>>> encoded;
    encoded.reserve(num_keys);
    for (int i = 0; i < num_keys; ++i) {
      encoded.push_back(std::make_pair(encode(vs[i]), vs[i]));
    }
    std::sort(encoded.begin(), encoded.end(),
              [](const auto& a, const auto& b) {
                if (a.first != b.first) return a.first < b.first;
                return a.second < b.second;
              });

    std::vector<std::pair<ZKey, uint64_t>> tree_entries;
    tree_entries.reserve(encoded.size());
    for (const auto& entry : encoded) {
      if (payloads_.empty() || entry.first != tree_entries.back().first) {
        payloads_.push_back({entry.second});
        tree_entries.push_back(
            std::make_pair(entry.first,
                           static_cast<uint64_t>(payloads_.size() - 1)));
      } else {
        payloads_.back().push_back(entry.second);
      }
    }

    if (!tree_entries.empty()) {
      btree_.bulk_load(tree_entries.data(),
                       static_cast<int>(tree_entries.size()));
    }
  }

  void insert(const Point<T>& key) {
    const ZKey z = encode(key);
    if (btree_.exists(z)) {
      auto payload_index = btree_.at(z);
      payloads_[payload_index].push_back(key);
      return;
    }
    payloads_.push_back({key});
    const uint64_t payload_index = static_cast<uint64_t>(payloads_.size() - 1);
    btree_.insert(z, payload_index);
  }

  bool exists(const Point<T>& key) {
    Point<T> found;
    return lookup(key, found);
  }

  bool lookup(const Point<T>& key, Point<T>& result) {
    const ZKey z = encode(key);
    if (!btree_.exists(z)) return false;

    auto payload_index = btree_.at(z);
    for (const auto& candidate : payloads_[payload_index]) {
      if (candidate == key) {
        result = candidate;
        return true;
      }
    }
    return false;
  }

  void rangeQuery(const Point<T>& min_key, const Point<T>& max_key,
                  std::vector<Point<T>>& result) {
    if (payloads_.empty()) return;
    const ZKey min_z = encode(min_key);
    const ZKey max_z = encode(max_key);
    const ZKey lower = std::min(min_z, max_z);
    const ZKey upper = std::max(min_z, max_z);
    auto z_range = btree_.rangeQuery(lower, upper);
    for (const auto& entry : z_range) {
      const auto& bucket = payloads_[entry.second];
      for (const auto& candidate : bucket) {
        if (in_range(candidate, min_key, max_key)) {
          result.push_back(candidate);
        }
      }
    }
  }

  size_t index_size() const {
    size_t total = 0;
    total += payloads_.size() * sizeof(std::vector<Point<T>>);
    for (const auto& bucket : payloads_) {
      total += bucket.capacity() * sizeof(Point<T>);
    }
    total += btree_.index_size();
    return total;
  }

 private:
 using ZKey = uint64_t;
  static constexpr uint32_t kPrecision = 4294967295u;  // 2^32 - 1
  struct Bounds {
    double min_x = 0.0;
    double min_y = 0.0;
    double span_x = 1.0;
    double span_y = 1.0;
  };

  void reset() {
    destroy_tree(btree_.root);
    btree_.root = new BTreeLeaf<ZKey, uint64_t>();
    payloads_.clear();
  }

  void destroy_tree(NodeBase* node) {
    if (!node) return;
    if (node->type == PageType::BTreeInner) {
      auto inner = static_cast<BTreeInner<ZKey>*>(node);
      for (int i = 0; i <= inner->count; ++i) {
        destroy_tree(inner->children[i]);
      }
    }
    delete node;
  }

  inline ZKey encode(const Point<T>& key) const {
    // If T is double, need to transform it to int first.
    if constexpr (std::is_same<T, double>::value) {
      const double norm_x =
          std::clamp((key.x - bounds_.min_x) / bounds_.span_x, 0.0, 1.0);
      const double norm_y =
          std::clamp((key.y - bounds_.min_y) / bounds_.span_y, 0.0, 1.0);
      return morton2D_64_encode(
          static_cast<uint32_t>(norm_x * kPrecision),
          static_cast<uint32_t>(norm_y * kPrecision));
    }
    return morton2D_64_encode(key.x, key.y);
  }

  inline bool in_range(const Point<T>& point, const Point<T>& min_key,
                       const Point<T>& max_key) const {
    return point.x >= min_key.x && point.x <= max_key.x &&
           point.y >= min_key.y && point.y <= max_key.y;
  }

  BTree<ZKey, uint64_t> btree_;
  std::vector<std::vector<Point<T>>> payloads_;
  Bounds bounds_;
};

#endif  // BTREE_ZORDER_HPP
