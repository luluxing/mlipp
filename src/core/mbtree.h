#ifndef __MBTREE_H__
#define __MBTREE_H__

#include "point.h"

// Store the 2D points in the way of MLIPP_kd
template <typename T>
class MBTree {
 public:
  MBTree() {}
  ~MBTree() {}

  void insert(const Point<T>& key) {

  }

  bool exists(const Point<T>& key) const {

  }

  void rangeQuery(const Point<T>& min_key,
                  const Point<T>& max_key,
                  std::vector<Point<T>>& result) {

  }

  int range_query(const Point<T>& lower, const Point<T>& upper, Point<T>* results) const {

  }

  int knn_query(const Point<T>& key, int k, Point<T>* results, double* distances) {

  }

  void bulk_load(Point<T>* vs, int num_keys) {

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
};

#endif // __MBTREE_H__