#ifndef __MBTREE_H__
#define __MBTREE_H__

#include <algorithm>
#include <cmath>
#include <sstream>
#include "point.h"

static const uint64_t PAGESIZE = 4*1024;

struct BNode {
  bool is_leaf;
  bool sort_by_x;
  uint16_t count;
};

template <typename T>
struct LeafNode : public BNode {
  static const uint16_t max_count = (PAGESIZE - sizeof(BNode)) / sizeof(Point<T>);
  Point<T> points[max_count];

  LeafNode(bool sort_by_x) {
    is_leaf = true;
    this->sort_by_x = sort_by_x;
    count = 0;
  }

  void load(Point<T>* vs, int begin, int num_keys) {
    for (int i = 0; i < num_keys; i++) {
      points[i] = vs[begin + i];
    }
    if (sort_by_x) {
      qsort(points, num_keys, sizeof(Point<T>), &compare_x<T>);
    } else {
      qsort(points, num_keys, sizeof(Point<T>), &compare_y<T>);
    }
    count = num_keys;
  }

  bool exists(const Point<T>& key) const {
    int lower = 0;
    int upper = count;
    do {
      int mid = ((upper - lower) / 2) + lower;
      if (get_dim(key, sort_by_x) < get_dim(points[mid], sort_by_x)) {
        upper = mid;
      } else if (get_dim(key, sort_by_x) > get_dim(points[mid], sort_by_x)) {
        lower = mid+1;
      } else {
        return PT_EQ(key, points[mid]);
      }
    } while (lower < upper);
    return PT_EQ(key, points[lower]);
  }
};

template <typename T>
struct InnerNode : public BNode {
  static const uint16_t max_count = (PAGESIZE - sizeof(BNode)) / (sizeof(T) +sizeof(BNode*));
  T keys[max_count];
  BNode* children[max_count];

  InnerNode(bool sort_by_x) {
    this->sort_by_x = sort_by_x;
    is_leaf = false;
    count = 0;
  }

  BNode* find_child(const Point<T>& key) const {
    int lower = 1;
    int upper = count;
    do {
      int mid = ((upper - lower) / 2) + lower;
      if (get_dim(key, sort_by_x) < keys[mid]) {
        upper = mid;
      } else if (get_dim(key, sort_by_x) > keys[mid]) {
        lower = mid+1;
      } else {
        return children[mid];
      }
    } while (lower < upper);
    return children[lower - 1];
  }
};

// Store the 2D points in the way of MLIPP_kd
template <typename T>
class MBTree {
 public:
  MBTree() { root_ = nullptr; }
  ~MBTree() {}

  void insert(const Point<T>& key) {

  }

  bool exists(const Point<T>& key) const {
    BNode* node = root_;
    while (true) {
      if (node->is_leaf) {
        return static_cast<LeafNode<T>*>(node)->exists(key);
      } else {
        InnerNode<T>* inner_node = static_cast<InnerNode<T>*>(node);
        node = inner_node->find_child(key);
      }
    }
    return false;

  }

  void rangeQuery(const Point<T>& min_key,
                  const Point<T>& max_key,
                  std::vector<Point<T>>& result) {

  }

  int range_query(const Point<T>& lower, const Point<T>& upper, Point<T>* results) const {

  }

  int knn_query(const Point<T>& key, int k, Point<T>* results, double* distances) {

  }

  // Points are not initially sorted
  // Duplicate x or y points are not allowed
  void bulk_load(Point<T>* vs, int num_keys) {
    root_ = bulk_load_helper(vs, 0, num_keys, true);

  }

  void show() const {

  }

  int max_depth() const {

  }

  void print_depth(int *max_depth_ret, double *avg_depth) const {

  }

  void print_stats() const {

  }

  size_t index_size(bool total=false, bool ignore_child=true) const {
    return 0;
  }

 private:
  BNode* root_;

  BNode* bulk_load_helper(Point<T>* vs, int begin, int num_keys, bool sort_by_x) {
    if (LeafNode<T>::max_count >= num_keys) {
      auto leaf_node = new LeafNode<T>(sort_by_x);
      leaf_node->load(vs, begin, num_keys);
      return leaf_node;
    }
    auto inner_node = new InnerNode<T>(sort_by_x);
    Point<T>* keys = vs + begin;
    if (sort_by_x)
      qsort(keys, num_keys, sizeof(Point<T>), &compare_x<T>);
    else
      qsort(keys, num_keys, sizeof(Point<T>), &compare_y<T>);
    
    int fanout = node_fanout(num_keys);
    int partition_size = std::ceil(1.0 * num_keys / fanout);
    int remaining = num_keys;
    for (int i = 0; i < fanout; i++) {
      if (sort_by_x)
        inner_node->keys[i] = keys[i * partition_size].x;
      else
        inner_node->keys[i] = keys[i * partition_size].y;
      auto child = bulk_load_helper(keys, i * partition_size,
                                    std::min(remaining, partition_size), !sort_by_x);
      inner_node->children[i] = child;
      inner_node->count++;
      remaining -= partition_size;
    }
    return inner_node;
  }

  int node_fanout(int num_keys) const {
    int node_num = num_keys % LeafNode<T>::max_count == 0 ?
                    num_keys / LeafNode<T>::max_count :
                    num_keys / LeafNode<T>::max_count + 1;
    int fanout = node_num;
    while (node_num > 1) {
      fanout = node_num;
      double r = 1.0 * node_num / InnerNode<T>::max_count;
      node_num = std::ceil(r);
    }
    return fanout;
  }
};

#endif // __MBTREE_H__