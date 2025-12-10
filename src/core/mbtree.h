#ifndef __MBTREE_H__
#define __MBTREE_H__

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stack>
#include "point.h"

namespace mbtree_internal {

static const uint64_t PAGESIZE = 4*1024;

template <typename T>
bool greater_than(T key, T query) {
  if (key == 0) return false;
  return key >= query;
}

template <typename T>
bool less_than(T key, T query) {
  if (key == 0) return false;
  return key <= query;
}

template <typename T>
bool overlap(T low, T up, T query_low, T query_up) {
  assert(!(low == 0 && up == 0));
  if (low == 0) return query_low <= up;
  if (up == 0) return query_up >= low;
  return !(up < query_low || query_up < low);
}

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

  void change_sort_axis() {
    sort_by_x = !sort_by_x;
    std::sort(points, points + count, [this](const Point<T>& a, const Point<T>& b) {
      return get_dim(a, sort_by_x) < get_dim(b, sort_by_x);
    });
  }

  bool is_full() const {
    return count == max_count;
  }

  LeafNode* split(const Point<T>& key, T& pivot) {
    auto new_leaf_node = new LeafNode<T>(sort_by_x);
    new_leaf_node->count = count - count / 2;
    count = count - new_leaf_node->count;
    for (int i = 0; i < new_leaf_node->count; i++) {
      new_leaf_node->points[i] = points[i + count];
    }
    pivot = get_dim(points[count - 1], sort_by_x);
    if (get_dim(key, sort_by_x) < pivot) {
      insert(key);
    } else {
      new_leaf_node->insert(key);
    }
    return new_leaf_node;
  }

  void insert(const Point<T>& key) {
    if (count == 0) {
      points[0] = key;
      count++;
      return;
    }
    int pos = lowerBound(get_dim(key, sort_by_x));
    memmove(points+pos+1, points+pos, sizeof(Point<T>)*(count-pos));
    points[pos] = key;
    count++;
  }

  int lowerBound(T k) {
    int lower=0;
    int upper=count;
    do {
       int mid = ((upper-lower)/2) + lower;
       if (k < get_dim(points[mid], sort_by_x)) {
          upper = mid;
       } else if (k > get_dim(points[mid], sort_by_x)) {
          lower = mid+1;
       } else {
          return mid;
       }
    } while (lower < upper);
    return lower;
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

  void get_data_within(const Point<T>& min_key, const Point<T>& max_key,
                       std::vector<Point<T>>& leaf_data) {
    for (int i = 0; i < count; i++) {
      if (points[i].x >= min_key.x && points[i].x <= max_key.x &&
          points[i].y >= min_key.y && points[i].y <= max_key.y) {
        leaf_data.push_back(points[i]);
      }
    }
  }

  void get_all_data(std::vector<Point<T>>& leaf_data) const {
    for (int i = 0; i < count; i++) {
      leaf_data.push_back(points[i]);
    }
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

  bool is_full() const {
    return count == max_count - 1;
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

  void replace_child(BNode* old_child, BNode* new_child) {
    for (int i = 0; i < count; i++) {
      if (children[i] == old_child) {
        children[i] = new_child;
        return;
      }
    }
  }

  int lower_bound(T pivot) {
    int lower = 0;
    int upper = count;
    do {
      int mid = ((upper - lower) / 2) + lower;
      if (pivot < keys[mid]) {
        upper = mid;
      } else if (pivot > keys[mid]) {
        lower = mid + 1;
      } else {
        return mid;
      }
    } while (lower < upper);
    return lower;
  }

  void insert(BNode* child, T pivot) {
    int pos = lower_bound(pivot);
    memmove(keys + pos + 1, keys + pos, sizeof(T) * (count-pos+1));
    memmove(children+pos+1, children+pos, sizeof(BNode*)*(count-pos+1));
    keys[pos] = pivot;
    children[pos] = child;
    std::swap(children[pos], children[pos+1]);
    count++;
  }

  bool child_within(int pos, const Point<T>& min_key, const Point<T>& max_key,
                    T& lower_x, T& upper_x, T& lower_y, T& upper_y) {
    T lower_key, upper_key;
    if (pos == 0) {
      lower_key = 0;
      upper_key = keys[0];
    } else if (pos == count) {
      lower_key = keys[count - 1];
      upper_key = 0;
    } else {
      lower_key = keys[pos - 1];
      upper_key = keys[pos];
    }
    if (sort_by_x) {
      lower_x = lower_key;
      upper_x = upper_key;
      return overlap(lower_x, upper_x, min_key.x, max_key.x);
    } else {
      lower_y = lower_key;
      upper_y = upper_key;
      return overlap(lower_y, upper_y, min_key.y, max_key.y);
    }
  }
};

// Store the 2D points in the way of MLIPP_kd
template <typename T>
class MBTree {
 public:
  MBTree() { root_ = new LeafNode<T>(true); }
  ~MBTree() {}

  // Since we change split axis alternatively, the tree cannot
  // be built from bottom to top like B+-tree, but can only be
  // built like kd-tree from top to bottom as we cannot group
  // entries in one leaf node initially.
  
  void insert(const Point<T>& key) {
    BNode* node = root_;
    BNode* parent = nullptr;
    while (!node->is_leaf) {
      InnerNode<T>* inner_node = static_cast<InnerNode<T>*>(node);
      int pos = inner_node->lower_bound(get_dim(key, inner_node->sort_by_x));
      node = inner_node->children[pos];
      parent = inner_node;
    }

    LeafNode<T>* leaf_node = static_cast<LeafNode<T>*>(node);
    if (leaf_node->is_full()) {
      T pivot;
      auto new_leaf_node = leaf_node->split(key, pivot);
      if (parent == nullptr) {
        root_ = make_parent(leaf_node, new_leaf_node, pivot);
      } else {
        InnerNode<T>* inner_node = static_cast<InnerNode<T>*>(parent);
        if (!inner_node->is_full()) {
          inner_node->insert(new_leaf_node, pivot);
        } else {
          auto new_inner_node = make_parent(leaf_node, new_leaf_node, pivot);
          inner_node->replace_child(leaf_node, new_inner_node);
          leaf_node->change_sort_axis();
          new_leaf_node->change_sort_axis();
        }
      }
    } else {
      leaf_node->insert(key);
    }

  }

  bool exists(const Point<T>& key) const {
    BNode* node = root_;
    while (true) {
      if (node->is_leaf) {
        return static_cast<LeafNode<T>*>(node)->exists(key);
      } else {
        InnerNode<T>* inner_node = static_cast<InnerNode<T>*>(node);
        int pos = inner_node->lower_bound(get_dim(key, inner_node->sort_by_x));
        node = inner_node->children[pos];
      }
    }
    return false;

  }

  void rangeQuery(const Point<T>& min_key,
                  const Point<T>& max_key,
                  std::vector<Point<T>>& result) {

    
    std::stack<NodeRange> s;
    s.push((NodeRange){0, 0, 0, 0, root_});
    while (!s.empty()) {
      auto node_range = s.top();
      s.pop();
      if (within_range(node_range, min_key, max_key)) {
        // Get all the leaf nodes in the range and add to result
        std::vector<Point<T>> leaf_data;
        get_subtree_data(node_range.node, leaf_data);
        // Append leaf_data to result
        result.insert(result.end(), leaf_data.begin(), leaf_data.end());
      } else {
        if (node_range.node->is_leaf) {
          // Find the data in the range and add to result
          std::vector<Point<T>> leaf_data;
          LeafNode<T>* leaf_node = static_cast<LeafNode<T>*>(node_range.node);
          leaf_node->get_data_within(min_key, max_key, leaf_data);
          // Append leaf_data to result
          result.insert(result.end(), leaf_data.begin(), leaf_data.end());
        } else {
          // Find the child nodes in the range and add to stack
          std::vector<NodeRange> child_nodes;
          get_children_within(node_range, min_key, max_key, child_nodes);
          for (auto child_node : child_nodes) {
            s.push(child_node);
          }
        }
      }
    }
  }

  int range_query(const Point<T>& lower, const Point<T>& upper, Point<T>* results) const {
    return 0;
  }

  int knn_query(const Point<T>& key, int k, Point<T>* results, double* distances) {
    return 0;
  }

  // Points are not initially sorted
  // Duplicate x or y points are not allowed
  void bulk_load(Point<T>* vs, int num_keys) {
    root_ = bulk_load_helper(vs, 0, num_keys, true);

  }

  // void show() const {

  // }

  // int max_depth() const {

  // }

  // void print_depth(int *max_depth_ret, double *avg_depth) const {

  // }

  // void print_stats() const {

  // }

  size_t index_size() const {
    std::stack<BNode*> s;
    s.push(root_);
    size_t size = 0;
    while (!s.empty()) {
      BNode* node = s.top();
      s.pop();
      size += sizeof(*node);
      if (!node->is_leaf) {
        InnerNode<T>* inner_node = static_cast<InnerNode<T>*>(node);
        for (int i = 0; i <= inner_node->count; i++) {
          s.push(inner_node->children[i]);
        }
      }
    }
    return size;
  }

 private:
  typedef struct {
    // Range of the node
    T lower_x;
    T upper_x;
    T lower_y;
    T upper_y;
    BNode* node;
  } NodeRange;
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
    for (int i = 0; i < fanout - 1; i++) {
      if (sort_by_x)
        inner_node->keys[i] = keys[(i + 1) * partition_size - 1].x;
      else
        inner_node->keys[i] = keys[(i + 1) * partition_size - 1].y;
      auto child = bulk_load_helper(keys, i * partition_size,
                                    std::min(remaining, partition_size), !sort_by_x);
      inner_node->children[i] = child;
      inner_node->count++;
      remaining -= partition_size;
      if (i == fanout - 2) {
        inner_node->children[i+1] = bulk_load_helper(keys, (i + 1) * partition_size,
                                    remaining, !sort_by_x);
      }
    }
    return inner_node;
  }

  BNode* make_parent(BNode* left, BNode* right, T pivot) {
    auto parent = new InnerNode<T>(left->sort_by_x);
    parent->keys[0] = pivot;
    parent->children[0] = left;
    parent->children[1] = right;
    parent->count = 1;
    return parent;
  }

  // TODO: check the fill ratio of each node after building the tree
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

  bool within_range(const NodeRange& node_range,
                    const Point<T>& min_key, const Point<T>& max_key) {
    
    return greater_than(node_range.lower_x, min_key.x) &&
           less_than(node_range.upper_x, max_key.x) &&
           greater_than(node_range.lower_y, min_key.y) &&
           less_than(node_range.upper_y, max_key.y);
  }

  void get_subtree_data(const BNode* node, std::vector<Point<T>>& leaf_data) {
    if (node->is_leaf) {
      const LeafNode<T>* leaf_node = static_cast<const LeafNode<T>*>(node);
      leaf_node->get_all_data(leaf_data);
    } else {
      const InnerNode<T>* inner_node = static_cast<const InnerNode<T>*>(node);
      for (int i = 0; i <= inner_node->count; i++) {
        get_subtree_data(inner_node->children[i], leaf_data);
      }
    }
  }

  void get_children_within(const NodeRange& node_range,
                           const Point<T>& min_key, const Point<T>& max_key,
                           std::vector<NodeRange>& child_nodes) {
    if (node_range.node->is_leaf) {
      return;
    }
    InnerNode<T>* inner_node = static_cast<InnerNode<T>*>(node_range.node);
    for (int i = 0; i <= inner_node->count; i++) {
      T lower_x = node_range.lower_x;
      T upper_x = node_range.upper_x;
      T lower_y = node_range.lower_y;
      T upper_y = node_range.upper_y;
      if (inner_node->child_within(i, min_key, max_key, lower_x, upper_x, lower_y, upper_y)) {
        child_nodes.push_back(NodeRange{lower_x, upper_x, lower_y, upper_y, inner_node->children[i]});
      }
    }
  }
};

} // namespace mbtree_internal

#endif // __MBTREE_H__