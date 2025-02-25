#ifndef __TWO_LIPP_H__
#define __TWO_LIPP_H__

#include <vector>
#include "lipp.h"
#include "point.h"

template<class T, bool USE_FMCD = true>
class LIPP_XY {
 public:
  LIPP_XY() {}

  ~LIPP_XY() {
    points.clear();
  }

  void bulk_load(Point<T>* vs, int num_keys) {}
  int knn_query(const Point<T>& key, int k, Point<T>* results, double* distances) {}
  int range_query(const Point<T>& lower, const Point<T>& upper, Point<T>* results) const {}
  void rangeQuery(const Point<T>& min_key,
    const Point<T>& max_key,
    std::vector<Point<T>>& result) {}
  void insert(const Point<T>& key) {}
  bool exists(const Point<T>& key) {}

 private:
  LIPP<T, uint64_t, USE_FMCD> lipp_x;
  LIPP<T, uint64_t, USE_FMCD> lipp_y;

  std::vector<Point<T>> points;
};

#endif // __TWO_LIPP_H__