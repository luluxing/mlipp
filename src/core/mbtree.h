#ifndef __MBTREE_H__
#define __MBTREE_H__

#include <algorithm>
#include "point.h"

static const uint64_t PAGESIZE = 4*1024;

struct Node {
  bool is_leaf;
  bool sort_by_x;
  uint16_t count;
};

template <typename T>
struct LeafNode : public Node {
  static const uint16_t max_count = (PAGESIZE - sizeof(Node)) / sizeof(Point<T>);
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

  T get_dim(const Point<T>& key) const {
    return sort_by_x ? key.x : key.y;
  }

  bool exists(const Point<T>& key) const {
    int lower = 0;
    int upper = count;
    do {
      int mid = ((upper - lower) / 2) + lower;
      if (get_dim(key) < get_dim(points[mid])) {
        upper = mid;
      } else if (get_dim(key) > get_dim(points[mid])) {
        lower = mid+1;
      } else {
        return true;
      }
    } while (lower < upper);
    return false;
  }
};

template <typename T>
struct InnerNode : public Node {
  static const uint16_t max_count = (PAGESIZE - sizeof(Node)) / (sizeof(T) +sizeof(Node*));
  T keys[max_count];
  Node* children[max_count];

  InnerNode(bool sort_by_x) {
    this->sort_by_x = sort_by_x;
    is_leaf = false;
    count = 0;
  }

  Node* find_child(const Point<T>& key) const {
    return children[0];
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
    Node* node = root_;
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
  void bulk_load(Point<T>* vs, int num_keys) {
    if (LeafNode<T>::max_count >= num_keys) {
      root_ = new LeafNode<T>(true);
      auto leaf_node = static_cast<LeafNode<T>*>(root_);
      leaf_node->load(vs, 0, num_keys);
    } else {
      root_ = new InnerNode<T>(true);
      
    }

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
  Node* root_;
};

#endif // __MBTREE_H__