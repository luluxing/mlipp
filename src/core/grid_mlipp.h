#ifndef GRID_MLIPP_H
#define GRID_MLIPP_H

#include <cassert>
#include <map>
#include <vector>
#include "mlipp_kd.h"

enum Partition {
  SPACE = 0,
  DATA = 1
};

// Always use space partition for simplicity
template<class T, bool USE_FMCD = true>
class GridMlipp {
 public:
  GridMlipp(int partition_x, int partition_y, Partition base)
  : partition_x_(partition_x), partition_y_(partition_y), base_(base) {
    if (partition_x <= 0 || partition_y <= 0 || base_ != Partition::SPACE) {
      fprintf(stderr, "Partition out of bound or base is not space.\n");
      exit(EXIT_FAILURE);
    }
    // fprintf(stdout, "Partition: %d, %d\n", partition_x_, partition_y_);
  }
  ~GridMlipp() {}

  void show() {
    for (const auto& pair : mlipp_map_) {
      pair.second->show();
    }
  }

  void bulk_load(Point<T>* vs, int num_keys) {
    // Partition the data or space based on the partition type
    if (base_ == Partition::DATA) {
      // Not implemented yet
    } else {
      assert(base_ == Partition::SPACE);
      // Obtain the minimum and maximum values for the axis of the given points
      T min_x = std::numeric_limits<T>::max();
      T max_x = std::numeric_limits<T>::min();
      T min_y = std::numeric_limits<T>::max();
      T max_y = std::numeric_limits<T>::min();
      for (int i = 0; i < num_keys; ++i) {
        min_x = std::min(min_x, vs[i].x);
        max_x = std::max(max_x, vs[i].x);
        min_y = std::min(min_y, vs[i].y);
        max_y = std::max(max_y, vs[i].y);
      }
      // Create partitions based on the axis
      step_x_ = (max_x - min_x) / partition_x_;
      step_y_ = (max_y - min_y) / partition_y_;
      last_x_ = min_x + (partition_x_ - 1) * step_x_;
      last_y_ = min_y + (partition_y_ - 1) * step_y_;
      for (int i = 0; i < partition_x_; ++i) {
        for (int j = 0; j < partition_y_; ++j) {
          T lower_bound_x = min_x + i * step_x_;
          T upper_bound_x = (i == partition_x_ - 1) ? std::numeric_limits<T>::max() : min_x + (i + 1) * step_x_;
          T lower_bound_y = min_y + j * step_y_;
          T upper_bound_y = (j == partition_y_ - 1) ? std::numeric_limits<T>::max() : min_y + (j + 1) * step_y_;
          Point<T> lower_left{lower_bound_x, lower_bound_y};
          Point<T> upper_right{upper_bound_x, upper_bound_y};
          mlipp_map_[lower_left] = new MLIPP_KD<T, USE_FMCD>(0.0, true);
          std::vector<Point<T>> partition_points;
          for (int j = 0; j < num_keys; ++j) {
            if ((vs[j].x >= lower_bound_x && vs[j].x < upper_bound_x) &&
                (vs[j].y >= lower_bound_y && vs[j].y < upper_bound_y)) {
              partition_points.push_back(vs[j]);
            }
          }
          mlipp_map_[lower_left]->bulk_load(partition_points.data(), partition_points.size());
        }
      }
    }
  }

  int knn_query(const Point<T>& key, int k, Point<T>* results, double* distances) {}

  int range_query(const Point<T>& lower, const Point<T>& upper, Point<T>* results) const {
    std::vector<Point<T>> temp_results;
    rangeQuery(lower, upper, temp_results);
    // Copy results to the output array
    int count = 0;
    for (const auto& pt : temp_results) {
      results[count++] = pt;
    }
    return count;
  }

  void rangeQuery(const Point<T>& min_key,
                  const Point<T>& max_key,
                  std::vector<Point<T>>& result) {
    // Break the range query into multiple partition queries and merge results
    for (auto it = mlipp_map_.begin(); it != mlipp_map_.end(); ++it) {
      auto mbr_lower_left = it->first;
      auto mbr_upper_right = Point<T>{mbr_lower_left.x + step_x_,
                                      mbr_lower_left.y + step_y_};
      if (mbr_lower_left.x > max_key.x || mbr_lower_left.y > max_key.y) {
        continue;
      } else if (mbr_upper_right.x <= min_key.x || mbr_upper_right.y <= min_key.y) {
        if (mbr_lower_left.x == last_x_ || mbr_lower_left.y == last_y_) {
          // Do nothing
        } else {
          continue;
        }
      }
      // Query the partition if it intersects with the range query
      std::vector<Point<T>> partition_results;
      partition_results.clear();
      it->second->rangeQuery(min_key, max_key, partition_results);
      result.insert(result.end(), partition_results.begin(), partition_results.end());
    }
  }

  void insert(const Point<T>& key) {}

  bool exists(const Point<T>& key) {
    for (auto it = mlipp_map_.begin(); it != mlipp_map_.end(); ++it) {
      auto mbr_lower_left = it->first;
      auto mbr_upper_right = Point<T>{mbr_lower_left.x + step_x_, mbr_lower_left.y + step_y_};
      if (mbr_lower_left.x > key.x || mbr_lower_left.y > key.y) {
        continue;
      } else if (mbr_upper_right.x <= key.x || mbr_upper_right.y <= key.y) {
        if (mbr_lower_left.x == last_x_ || mbr_lower_left.y == last_y_) {
            // Do nothing
        } else {
          continue;
        }
      }
      if (!it->second->exists(key)) continue;
      else return true;
    }
    return false;
  }

  size_t index_size() const {
    size_t total_size = 0;
    for (const auto& pair : mlipp_map_) {
      total_size += pair.second->index_size();
    }
    return total_size;
  }

 private:
  int partition_x_;
  int partition_y_;
  T last_x_;
  T last_y_;
  T step_x_;
  T step_y_;
  Partition base_ = Partition::SPACE;
  std::map<Point<T>, MLIPP_KD<T, USE_FMCD>*> mlipp_map_;
};

#endif // GRID_MLIPP_H