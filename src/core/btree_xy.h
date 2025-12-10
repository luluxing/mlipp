#ifndef BTREE_XY_HPP
#define BTREE_XY_HPP

#include "btree_olc.h"
#include "point.h"

using namespace btree_olc_internal;

template<class T>
class BTreeXY {
  public:
  BTreeXY() {}

  ~BTreeXY() {
    pts_sort_x.clear();
    pts_sort_y.clear();
  }

  void bulk_load(Point<T>* vs, int num_keys) {
    for (int i = 0; i < num_keys; ++i) {
      pts_sort_x.push_back(vs[i]);
      pts_sort_y.push_back(vs[i]);
    }
    // Sort the points by x and y
    std::sort(pts_sort_x.begin(), pts_sort_x.end(), 
        [](const auto& a, const auto& b) {
          return a.x < b.x;
        });
    std::sort(pts_sort_y.begin(), pts_sort_y.end(),
        [](const auto& a, const auto& b) {
          return a.y < b.y;
        });

    // Populate btree_x
    std::vector<std::pair<T, uint64_t>> x_coords;
    int x_num = 0;
    for (int i = 0; i < num_keys; ++i) {
      if (i == 0 || pts_sort_x[i].x != pts_sort_x[i - 1].x) {
        x_coords.push_back(std::make_pair(pts_sort_x[i].x, i));
        x_num++;
      }
    }
    btree_x.bulk_load(x_coords.data(), x_num);
    // Populate btree_y
    std::vector<std::pair<T, uint64_t>> y_coords;
    int y_num = 0;
    for (int i = 0; i < num_keys; ++i) {
      if (i == 0 || pts_sort_y[i].y != pts_sort_y[i - 1].y) {
        y_coords.push_back(std::make_pair(pts_sort_y[i].y, i));
        y_num++;
      }
    }

    btree_y.bulk_load(y_coords.data(), y_num);
  }

  int knn_query(const Point<T>& key, int k, Point<T>* results, double* distances) {
    fprintf(stderr, "KNN query not supported in BTreeXY\n");
    return 0;
  }

  void rangeQuery(const Point<T>& min_key,
                  const Point<T>& max_key,
                  std::vector<Point<T>>& result) {
    // We only query the smaller range of x or y
    T x_span = max_key.x - min_key.x;
    T y_span = max_key.y - min_key.y;
    if (x_span < y_span) {
      // Query x
      std::vector<std::pair<T, uint64_t>> x_result = btree_x.rangeQuery(min_key.x, max_key.x);
      // Filter by y
      for (auto& x_ind : x_result) {
        auto& pt = pts_sort_x[x_ind.second];
        for (int xi = x_ind.second; xi < pts_sort_x.size(); ++xi) {
          auto temp = pts_sort_x[xi];
          if (pt.x != temp.x) break;
          if (pt.y >= min_key.y && pt.y <= max_key.y) {
            result.push_back(temp);
          } else if (pt.y > max_key.y) {
            break;
          }
        }
      }
    } else {
      // Query y
      std::vector<std::pair<T,uint64_t>> y_result = btree_y.rangeQuery(min_key.y, max_key.y);
      // Filter by x
      for (auto& y_ind : y_result) {
        auto& pt = pts_sort_y[y_ind.second];
        for (int yi = y_ind.second; yi < pts_sort_y.size(); ++yi) {
          auto temp = pts_sort_y[yi];
          if (pt.y != temp.y) break;
          if (pt.x >= min_key.x && pt.x <= max_key.x) {
            result.push_back(temp);
          } else if (pt.x > max_key.x) {
            break;
          }
        }
      }
    }
  }

  void insert(const Point<T>& key) {
    // Not supported
    fprintf(stderr, "Insert not supported in BTreeXY\n");

  }

  bool exists(const Point<T>& key) {
    if (!btree_x.exists(key.x)) return false;
    if (!btree_y.exists(key.y)) return false;

    auto x_pos = btree_x.at(key.x);
    for (; x_pos < pts_sort_x.size(); ++x_pos) {
      if (pts_sort_x[x_pos].x > key.x) break;
      if (pts_sort_x[x_pos].y == key.y) return true;
    }
    return false;
  }

  size_t index_size() const {
    auto s = btree_x.index_size() + btree_y.index_size();
    s += pts_sort_x.size() * sizeof(Point<T>);
    s += pts_sort_y.size() * sizeof(Point<T>);
    return s;
  }

 private:
  BTree<T, uint64_t> btree_x;
  BTree<T, uint64_t> btree_y;

  // In the case of same x or y coordinates,
  // we only store the location of the first one
  // in btree.
  std::vector<Point<T>> pts_sort_x;
  std::vector<Point<T>> pts_sort_y;
};

#endif