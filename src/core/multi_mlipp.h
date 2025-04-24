#ifndef MULTI_MLIPP_H
#define MULTI_MLIPP_H

#include <map>
#include <vector>
#include "mlipp_kd.h"

enum Axis {
  X_AXIS = 0,
  Y_AXIS = 1
};

enum Partition {
  SPACE = 0,
  DATA = 1
};

template<class T, bool USE_FMCD = true>
class MultiMlipp {
 public:
  MultiMlipp(int partition, Axis ax, Partition base)
  : partition_(partition), part_axis_(ax), base_(base) {
    if (partition <= 0 || partition > 5) {
      fprintf(stderr, "Partition out of bound.\n");
      exit(EXIT_FAILURE);
    }
  }
  ~MultiMlipp() {}

  void bulk_load(Point<T>* vs, int num_keys) {
    // Partition the data or space based on the partition type
    if (base_ == Partition::DATA) {
      std::vector<T> delimiters;
      // Sort the data given the axis into partitions
      if (part_axis_ == Axis::X_AXIS) {
        // Sort by X-axis
        std::sort(vs, vs + num_keys, [](const Point<T>& a, const Point<T>& b) {
          return a.x < b.x;
        });
        int remaining = num_keys;
        int num_keys_per_partition = num_keys / partition_;
        for (int i = 0 ; i < partition_; ++i) {
          delimiters.push_back(vs[i * num_keys_per_partition].x);
          mlipp_map_[delimiters.back()] = new MLIPP_KD<T, USE_FMCD>(0.0, true);
          if (i == partition_ - 1) {
            // Last partition may contain remaining elements
            mlipp_map_[delimiters.back()]->bulk_load(vs + i * num_keys_per_partition,
                                                     remaining);
          } else {
            mlipp_map_[delimiters.back()]->bulk_load(vs + i * num_keys_per_partition,
                                                     num_keys_per_partition);
          }
          remaining -= num_keys_per_partition;
        }
      } else {
        // Sort by Y-axis
        std::sort(vs, vs + num_keys, [](const Point<T>& a, const Point<T>& b) {
          return a.y < b.y;
        });
        int remaining = num_keys;
        int num_keys_per_partition = num_keys / partition_;
        for (int i = 0 ; i < partition_; ++i) {
          delimiters.push_back(vs[i * num_keys_per_partition].y);
          // print delimiter
          mlipp_map_[delimiters.back()] = new MLIPP_KD<T, USE_FMCD>(0.0, true);
          if (i == partition_ - 1) {
            // Last partition may contain remaining elements
            mlipp_map_[delimiters.back()]->bulk_load(vs + i * num_keys_per_partition,
                                                     remaining);
          } else {
            mlipp_map_[delimiters.back()]->bulk_load(vs + i * num_keys_per_partition,
                                                     num_keys_per_partition);
          }
          remaining -= num_keys_per_partition;
        }
      }
    } else {
      assert(base_ == Partition::SPACE);
      // Obtain the minimum and maximum values for the axis of the given points
      T min_val = std::numeric_limits<T>::max();
      T max_val = std::numeric_limits<T>::min();
      for (int i = 0; i < num_keys; ++i) {
        if (part_axis_ == Axis::X_AXIS) {
          min_val = std::min(min_val, vs[i].x);
          max_val = std::max(max_val, vs[i].x);
        } else {
          min_val = std::min(min_val, vs[i].y);
          max_val = std::max(max_val, vs[i].y);
        }
      }
      // Create partitions based on the axis
      T step = (max_val - min_val) / partition_;
      for (int i = 0; i < partition_; ++i) {
        T lower_bound = min_val + i * step;
        T upper_bound = (i == partition_ - 1) ? std::numeric_limits<T>::max() : min_val + (i + 1) * step;
        mlipp_map_[lower_bound] = new MLIPP_KD<T, USE_FMCD>(0.0, true);
        std::vector<Point<T>> partition_points;
        for (int j = 0; j < num_keys; ++j) {
          if ((part_axis_ == Axis::X_AXIS && vs[j].x >= lower_bound && vs[j].x < upper_bound) ||
              (part_axis_ == Axis::Y_AXIS && vs[j].y >= lower_bound && vs[j].y < upper_bound)) {
            partition_points.push_back(vs[j]);
          }
        }
        mlipp_map_[lower_bound]->bulk_load(partition_points.data(), partition_points.size());
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
    T min_val = part_axis_ == Axis::X_AXIS ? min_key.x : min_key.y;
    T max_val = part_axis_ == Axis::X_AXIS ? max_key.x : max_key.y;
    auto it_low = mlipp_map_.lower_bound(min_val);
    if (it_low != mlipp_map_.end() && it_low->first == min_val) {
      // Do nothing
    } else {
      if (it_low != mlipp_map_.begin()) {
        --it_low; // Move to the previous partition if the exact key is not found
      }
    }
    auto it_up = mlipp_map_.upper_bound(max_val);
    for (auto it = it_low; it != it_up; ++it) {
      // Query each partition
      std::vector<Point<T>> partition_results;
      it->second->rangeQuery(min_key, max_key, partition_results);
      result.insert(result.end(), partition_results.begin(), partition_results.end());
    }
  }

  void insert(const Point<T>& key) {}

  bool exists(const Point<T>& key) {
    T target_pt = part_axis_ == Axis::X_AXIS ? key.x : key.y;
    // Iterate over the leaf and find the partition that contains the key
    auto it = mlipp_map_.lower_bound(target_pt);
    if (it != mlipp_map_.end() && it->first == target_pt) {
      return it->second->exists(key);
    }
    if (it != mlipp_map_.begin()) {
      --it; // Move to the previous partition if the exact key is not found
      return it->second->exists(key);
    }
    return false; // Key not found in any partition
  }

  size_t index_size() const {}

 private:
  int partition_;
  Axis part_axis_;
  Partition base_;
  std::map<double, MLIPP_KD<T, USE_FMCD>*> mlipp_map_;
};

#endif // MULTI_MLIPP_H