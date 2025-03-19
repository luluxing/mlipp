#ifndef __MLIPP_Z_H__
#define __MLIPP_Z_H__

#include <stdint.h>
#include <math.h>
#include <limits>
#include <cstdio>
#include <stack>
#include <queue>
#include <vector>
#include <cstring>
#include <sstream>

#include "lipp_base.h"
#include "lipp_utils.h"
#include "morton.h"
#include "point.h"

using namespace libmorton;

uint32_t precision = 4294967295; // 2^32 - 1
// uint32_t precision = 8388607; // 2^23 - 1


template<class T, bool USE_FMCD = true>
class MLIPP_Z
{
  typedef uint_fast64_t Z;

  inline int compute_gap_count(int size) {
    if (size >= 1000000) return 1;
    if (size >= 100000) return 2;
    return 5;
  }

  inline Z encode(const Point<T>& key) const {
    // If T is double, need to transform it to int first.
    if constexpr (std::is_same<T, double>::value) {
      return morton2D_64_encode(
        static_cast<uint32_t>(key.x * precision),
        static_cast<uint32_t>(key.y * precision));
    }
    return morton2D_64_encode(key.x, key.y);
  }

  inline Point<T> decode(const Z& zkey) const {
    uint_fast32_t x_uint;
    uint_fast32_t y_uint;
    morton2D_64_decode(zkey, x_uint, y_uint);
    Point<T> key;
    key.x = static_cast<int>(x_uint);
    key.y = static_cast<int>(y_uint);
    return key;
  }

  struct Node;
  inline int PREDICT_POS(Node* node, Z zkey) const {
    double v = node->model.predict_double(zkey);
    if (v > std::numeric_limits<int>::max() / 2) {
      return node->num_items - 1;
    }
    if (v < 0) {
      return 0;
    }
    return std::min(node->num_items - 1, static_cast<int>(v));
  }

  static void remove_last_bit(bitmap_t& bitmap_item) {
    bitmap_item -= 1 << BITMAP_NEXT_1(bitmap_item);
  }

  const double BUILD_LR_REMAIN;
  const bool QUIET;

  struct {
    long long fmcd_success_times = 0;
    long long fmcd_broken_times = 0;
  } stats;

public:

  MLIPP_Z(double BUILD_LR_REMAIN = 0, bool QUIET = true)
    : BUILD_LR_REMAIN(BUILD_LR_REMAIN), QUIET(QUIET) {
    /*{
      std::vector<Node*> nodes;
      for (int _ = 0; _ < 1e7; _ ++) {
        Node* node = build_tree_two(Point(0, 0), Point(1, 1));
        nodes.push_back(node);
      }
      for (auto node : nodes) {
        destroy_tree(node);
      }
      if (!QUIET) {
        printf("initial memory pool size = %lu\n", pending_two.size());
      }
    }
    if (USE_FMCD && !QUIET) {
      printf("enable FMCD\n");
    }*/

    root = build_tree_none();
  }
  ~MLIPP_Z() {
    destroy_tree(root);
    root = NULL;
    destory_pending();
  }

  void insert(const Point<T>& key) {
    root = insert_tree(root, key);
  }
  bool exists(const Point<T>& key) {
    const Z zkey = encode(key);
    Node* node = root;
    while (true) {
      int pos = PREDICT_POS(node, zkey);
      if (BITMAP_GET(node->none_bitmap, pos) == 1) {
        return false;
      } else if (BITMAP_GET(node->child_bitmap, pos) == 0) {
        return PT_EQ(node->items[pos].comp.data, key);
      } else {
        node = node->items[pos].comp.child;
      }
    }
  }
  void rangeQuery(const Point<T>& min_key, const Point<T>& max_key, std::vector<Point<T>> &result) {
    const Z min_zkey = encode(min_key);
    const Z max_zkey = encode(max_key);
    rangeQueryInternal(root, min_key, max_key, min_zkey, max_zkey, result);
    return;
  }
  // Find the keys which are in range [lower, upper], returns the number of found keys.
  int range_query(const Point<T>& lower, const Point<T>& upper, Point<T>* results) const {
    const Z zlower = encode(lower);
    const Z zupper = encode(upper);
    return range_core<false, false>(results, 0, root, lower, upper, zlower, zupper);
  }
  void bulk_load(Point<T>* vs, int num_keys) {
    if (num_keys == 0) {
      destroy_tree(root);
      root = build_tree_none();
      return;
    }
    if (num_keys == 1) {
      destroy_tree(root);
      root = build_tree_none();
      insert(vs[0]);
      return;
    }
    if (num_keys == 2) {
      destroy_tree(root);
      root = build_tree_two(vs[0], encode(vs[0]), vs[1], encode(vs[1]));
      return;
    }

    RT_ASSERT(num_keys > 2);
    /*for (int i = 1; i < num_keys; i ++) {
      RT_ASSERT(vs[i].x > vs[i-1].x);
    }*/

    destroy_tree(root);
    root = build_tree_bulk(vs, num_keys, false);
  }

  void show() const {
    // printf("============= SHOW MLIPP_Z ================\n");
    fprintf(stdout, "Max depth: %d\n", max_depth());

    std::queue<Node*> s;
    s.push(root);
    while (!s.empty()) {
      // Node* node = s.top();
      Node* node = s.front();
      s.pop();

      // printf("Node(%p, a = %lf, b = %lf, num_items = %d)", node, node->model.a, node->model.b, node->num_items);
      // printf("[");
      int x = 1;
      for (int i = 0; i < node->num_items; i ++) {
        // if (!x) {
        //   printf(", ");
        // }
        x = 0;
        if (BITMAP_GET(node->none_bitmap, i) == 1) {
          // printf("None");
        } else if (BITMAP_GET(node->child_bitmap, i) == 0) {
          std::stringstream s;
          s << node->items[i].comp.data;
          printf("%s ", s.str().c_str());
        } else {
          // printf("Child(%p)", node->items[i].comp.child);
          s.push(node->items[i].comp.child);
        }
      }
      printf("\n");
    }
  }
  int max_depth() const {
    std::stack<Node*> s;
    s.push(root);

    int max_depth = 1;
    while (!s.empty()) {
      Node* node = s.top(); s.pop();
      for (int i = 0; i < node->num_items; i ++) {
        if (BITMAP_GET(node->child_bitmap, i) == 1) {
          s.push(node->items[i].comp.child);
        }
      }
    }
    return max_depth;
  }
  void print_depth(int *max_depth_ret, double *avg_depth) const {
    std::stack<Node*> s;
    std::stack<int> d;
    s.push(root);
    d.push(1);

    int max_depth = 1;
    int sum_depth = 0, sum_nodes = 0;
    while (!s.empty()) {
      Node* node = s.top(); s.pop();
      int depth = d.top(); d.pop();
      for (int i = 0; i < node->num_items; i ++) {
        if (BITMAP_GET(node->child_bitmap, i) == 1) {
          s.push(node->items[i].comp.child);
          d.push(depth + 1);
        } else if (BITMAP_GET(node->none_bitmap, i) != 1) {
          max_depth = std::max(max_depth, depth);
          sum_depth += depth;
          sum_nodes ++;
        }
      }
    }

    // printf("max_depth = %d, avg_depth = %.2lf\n", max_depth, double(sum_depth) / double(sum_nodes));
    *max_depth_ret = max_depth;
    *avg_depth = double(sum_depth) / double(sum_nodes);
    return;
  }
  void verify() const {
    std::stack<Node*> s;
    s.push(root);

    while (!s.empty()) {
      Node* node = s.top(); s.pop();
      int sum_size = 0;
      for (int i = 0; i < node->num_items; i ++) {
        if (BITMAP_GET(node->child_bitmap, i) == 1) {
          s.push(node->items[i].comp.child);
          sum_size += node->items[i].comp.child->size;
        } else if (BITMAP_GET(node->none_bitmap, i) != 1) {
          sum_size ++;
        }
      }
      RT_ASSERT(sum_size == node->size);
    }
  }
  void print_stats() const {
    printf("======== Stats ===========\n");
    if (USE_FMCD) {
      printf("\t fmcd_success_times = %lld\n", stats.fmcd_success_times);
      printf("\t fmcd_broken_times = %lld\n", stats.fmcd_broken_times);
    }
  }
  size_t index_size(bool total=false, bool ignore_child=true) const {
    std::stack<Node*> s;
    s.push(root);
  
    size_t size = 0;
    while (!s.empty()) {
      Node* node = s.top(); s.pop();
      bool has_child = false;
      if(ignore_child == false) {
        size += sizeof(*node);
      }
      for (int i = 0; i < node->num_items; i ++) {
        if (ignore_child == true) {
          size += sizeof(Item);
          has_child = true;
        } else {
          if (total) size += sizeof(Item);
        }
        if (BITMAP_GET(node->child_bitmap, i) == 1) {
          if (!total) size += sizeof(Item);
          s.push(node->items[i].comp.child);
        }
      }
      if (ignore_child == true && has_child) {
        size += sizeof(*node);
      }
    }
    return size;
  }

private:
  struct Node;
  struct Item
  {
    union {
      Point<T> data;
      Node* child;
    } comp;
  };
  struct Node
  {
    int is_two; // is special node for only two keys
    int build_size; // tree size (include sub nodes) when node created
    int size; // current tree size (include sub nodes)
    int fixed; // fixed node will not trigger rebuild
    int num_inserts, num_insert_to_data;
    int num_items; // size of items
    LinearModel<Z> model;
    Item* items;
    bitmap_t* none_bitmap; // 1 means None, 0 means Data or Child
    bitmap_t* child_bitmap; // 1 means Child. will always be 0 when none_bitmap is 1
  };

  Node* root;
  std::stack<Node*> pending_two;

  std::allocator<Node> node_allocator;
  Node* new_nodes(int n)
  {
    Node* p = node_allocator.allocate(n);
    RT_ASSERT(p != NULL && p != (Node*)(-1));
    return p;
  }
  void delete_nodes(Node* p, int n)
  {
    node_allocator.deallocate(p, n);
  }

  std::allocator<Item> item_allocator;
  Item* new_items(int n)
  {
    Item* p = item_allocator.allocate(n);
    RT_ASSERT(p != NULL && p != (Item*)(-1));
    return p;
  }
  void delete_items(Item* p, int n)
  {
    item_allocator.deallocate(p, n);
  }

  std::allocator<bitmap_t> bitmap_allocator;
  bitmap_t* new_bitmap(int n)
  {
    bitmap_t* p = bitmap_allocator.allocate(n);
    RT_ASSERT(p != NULL && p != (bitmap_t*)(-1));
    return p;
  }
  void delete_bitmap(bitmap_t* p, int n)
  {
    bitmap_allocator.deallocate(p, n);
  }

  // SATISFY_LOWER = true means all the keys in the subtree of `node` are no less than to `lower`.
  // SATISFY_UPPER = true means all the keys in the subtree of `node` are no greater than to `upper`.
  template<bool SATISFY_LOWER, bool SATISFY_UPPER>
  int range_core(Point<T>* results, int pos, Node* node, 
    const Point<T>& lower, const Point<T>& upper,
    const Z zlower, const Z zupper) const
  {
    if constexpr (SATISFY_LOWER && SATISFY_UPPER) {
      int bit_pos = 0;
      const bitmap_t* none_bitmap = node->none_bitmap;
      while (bit_pos < node->num_items) {
        bitmap_t not_none = ~(*none_bitmap);
        while (not_none) {
          int latest_pos = BITMAP_NEXT_1(not_none);
          not_none ^= 1 << latest_pos;

          int i = bit_pos + latest_pos;
          if (BITMAP_GET(node->child_bitmap, i) == 0) {
            if (node->items[i].comp.data.x >= lower.x 
             && node->items[i].comp.data.x <= upper.x
             && node->items[i].comp.data.y >= lower.y
             && node->items[i].comp.data.y <= upper.y)
            {
              results[pos] = node->items[i].comp.data;
              pos++;
            }
          } else {
            pos = range_core<true, true>(results, pos, node->items[i].comp.child, 
              lower, upper, zlower, zupper);
          }
        }

        bit_pos += BITMAP_WIDTH;
        none_bitmap ++;
      }
      return pos;
    } else {
      int lower_pos = SATISFY_LOWER ? -1 : PREDICT_POS(node, zlower);
      int upper_pos = SATISFY_UPPER ? node->num_items : PREDICT_POS(node, zupper);
      if constexpr (!SATISFY_LOWER) {
        if (BITMAP_GET(node->none_bitmap, lower_pos) == 0) {
          if (BITMAP_GET(node->child_bitmap, lower_pos) == 0) {
            do {
              if (encode(node->items[lower_pos].comp.data) < zlower) break;
              if constexpr (!SATISFY_UPPER) {
                if (encode(node->items[lower_pos].comp.data) > zupper) break;
              }
              if (node->items[lower_pos].comp.data.x >= lower.x 
               && node->items[lower_pos].comp.data.x <= upper.x
               && node->items[lower_pos].comp.data.y >= lower.y
               && node->items[lower_pos].comp.data.y <= upper.y)
              {
                results[pos] = node->items[lower_pos].comp.data;
                pos++;
              }
            } while (false);
          } else {
            if (lower_pos < upper_pos) {
              pos = range_core<false, true>(results, pos, node->items[lower_pos].comp.child, 
                lower, upper, zlower, zupper);
            } else {
              pos = range_core<false, false>(results, pos, node->items[lower_pos].comp.child, 
                lower, upper, zlower, zupper);
            }
          }
        }
      }
      {
        int bit_pos = (lower_pos + 1) / BITMAP_WIDTH * BITMAP_WIDTH;
        const bitmap_t* none_bitmap = node->none_bitmap + bit_pos / BITMAP_WIDTH;
        while (bit_pos < upper_pos) {
          bitmap_t not_none = ~(*none_bitmap);
          while (not_none) {
            int latest_pos = BITMAP_NEXT_1(not_none);
            not_none ^= 1 << latest_pos;

            int i = bit_pos + latest_pos;
            if (i <= lower_pos) continue;
            if (i >= upper_pos) break;

            if (BITMAP_GET(node->child_bitmap, i) == 0) {
              if (node->items[i].comp.data.x >= lower.x 
               && node->items[i].comp.data.x <= upper.x
               && node->items[i].comp.data.y >= lower.y
               && node->items[i].comp.data.y <= upper.y)
              {
                results[pos] = node->items[i].comp.data;
                pos++;
              }
            } else {
              pos = range_core<true, true>(results, pos, node->items[i].comp.child, 
                lower, upper, zlower, zupper);
            }
          }

          bit_pos += BITMAP_WIDTH;
          none_bitmap ++;
        }
      }
      if constexpr (!SATISFY_UPPER) {
        if (lower_pos < upper_pos) {
          if (BITMAP_GET(node->none_bitmap, upper_pos) == 0) {
            if (BITMAP_GET(node->child_bitmap, upper_pos) == 0) {
              if (encode(node->items[upper_pos].comp.data) <= zupper) {
                if (node->items[upper_pos].comp.data.x >= lower.x 
                 && node->items[upper_pos].comp.data.x <= upper.x
                 && node->items[upper_pos].comp.data.y >= lower.y
                 && node->items[upper_pos].comp.data.y <= upper.y)
                {
                  results[pos] = node->items[upper_pos].comp.data;
                  pos++;
                }
              }
            } else {
              pos = range_core<true, false>(results, pos, node->items[upper_pos].comp.child, 
                lower, upper, zlower, zupper);
            }
          }
        }
      }
      return pos;
    }
  }

  void rangeQueryInternal(Node* root,
    const Point<T>& min_key, const Point<T>& max_key,
    const Z& min_zkey, const Z& max_zkey,
    std::vector<Point<T>> &result) {
    
    int min_pos = PREDICT_POS(root, min_zkey);
    int max_pos = PREDICT_POS(root, max_zkey);
    for (int i = min_pos; i < max_pos + 1; i ++) {
      if (BITMAP_GET(root->none_bitmap, i) == 0) {
        if (BITMAP_GET(root->child_bitmap, i) == 0) {
          if (root->items[i].comp.data.x >= min_key.x 
           && root->items[i].comp.data.x <= max_key.x
           && root->items[i].comp.data.y >= min_key.y
           && root->items[i].comp.data.y <= max_key.y)
            result.push_back(root->items[i].comp.data);
        } else {
          rangeQueryInternal(root->items[i].comp.child,
                     min_key, max_key,
                     min_zkey, max_zkey,
                     result);
        }
      }
    }
    return;
  }

  /// build an empty tree
  Node* build_tree_none()
  {
    Node* node = new_nodes(1);
    node->is_two = 0;
    node->build_size = 0;
    node->size = 0;
    node->fixed = 0;
    node->num_inserts = node->num_insert_to_data = 0;
    node->num_items = 1;
    node->model.a = node->model.b = 0;
    node->items = new_items(1);
    node->none_bitmap = new_bitmap(1);
    node->none_bitmap[0] = BITMAP_ZERO;
    BITMAP_SET(node->none_bitmap, BITMAP_ZERO);
    node->child_bitmap = new_bitmap(1);
    node->child_bitmap[0] = BITMAP_ZERO;

    return node;
  }
  /// build a tree with two keys
  Node* build_tree_two(Point<T> key1, Z zkey1, Point<T> key2, Z zkey2)
  {
    if (zkey1 > zkey2) {
      std::swap(key1, key2);
      std::swap(zkey1, zkey2);
    }
    RT_ASSERT(zkey1 < zkey2);
    static_assert(BITMAP_WIDTH == 8);

    Node* node = NULL;
    if (pending_two.empty()) {
      node = new_nodes(1);
      node->is_two = 1;
      node->build_size = 2;
      node->size = 2;
      node->fixed = 0;
      node->num_inserts = node->num_insert_to_data = 0;

      node->num_items = 8;
      node->items = new_items(node->num_items);
      node->none_bitmap = new_bitmap(1);
      node->child_bitmap = new_bitmap(1);
      node->none_bitmap[0] = BITMAP_ONE;
      node->child_bitmap[0] = BITMAP_ZERO;
    } else {
      node = pending_two.top(); pending_two.pop();
    }

    const long double mid1_key = zkey1;
    const long double mid2_key = zkey2;

    const double mid1_target = node->num_items / 3;
    const double mid2_target = node->num_items * 2 / 3;

    node->model.a = (mid2_target - mid1_target) / (mid2_key - mid1_key);
    node->model.b = mid1_target - node->model.a * mid1_key;
    RT_ASSERT(isfinite(node->model.a));
    RT_ASSERT(isfinite(node->model.b));

    { // insert key1&value1
      int pos = PREDICT_POS(node, zkey1);
      RT_ASSERT(BITMAP_GET(node->none_bitmap, pos) == 1);
      BITMAP_CLEAR(node->none_bitmap, pos);
      node->items[pos].comp.data = key1;
    }
    { // insert key2&value2
      int pos = PREDICT_POS(node, zkey2);
      RT_ASSERT(BITMAP_GET(node->none_bitmap, pos) == 1);
      BITMAP_CLEAR(node->none_bitmap, pos);
      node->items[pos].comp.data = key2;
    }

    return node;
  }
  /// bulk build, _keys must be sorted in asc order.
  Node* build_tree_bulk(Point<T>* _keys, int _size, bool sorted)
  {
    if (USE_FMCD) {
      return build_tree_bulk_fmcd(_keys, _size, sorted);
    } else {
      return build_tree_bulk_fast(_keys, _size, sorted);
    }
  }
  /// bulk build, _keys must be sorted in asc order.
  /// split keys into three parts at each node.
  Node* build_tree_bulk_fast(Point<T>* _keys, int _size, bool sorted)
  {
    Z* _zkeys = new Z[_size];
    for (int i = 0; i < _size; ++i)
      _zkeys[i] = encode(_keys[i]);
    if (!sorted)
    {
      std::pair<Z, Point<T>>* pairs = new std::pair<Z, Point<T>>[_size];
      for (int i = 0; i < _size; ++i )
        pairs[i] = std::make_pair(_zkeys[i], _keys[i]);
      std::sort(pairs, pairs + _size, [](std::pair<Z, Point<T>> &left, std::pair<Z, Point<T>> &right) {
        return left.first < right.first;
      });
      for (int i = 0; i < _size; ++i )
      {
        _zkeys[ i ] = pairs[i].first;
        _keys[ i ] = pairs[i].second;
      }
      delete[] pairs;
    }

    RT_ASSERT(_size > 1);

    typedef struct {
      int begin;
      int end;
      int level; // top level = 1
      Node* node;
    } Segment;
    std::stack<Segment> s;

    Node* ret = new_nodes(1);
    s.push((Segment){0, _size, 1, ret});

    while (!s.empty()) {
      const int begin = s.top().begin;
      const int end = s.top().end;
      const int level = s.top().level;
      Node* node = s.top().node;
      s.pop();

      RT_ASSERT(end - begin >= 2);
      if (end - begin == 2) {
        Node* _ = build_tree_two(_keys[begin], _zkeys[begin], _keys[begin+1], _zkeys[begin+1]);
        memcpy(node, _, sizeof(Node));
        delete_nodes(_, 1);
      } else {
        Point<T>* keys = _keys + begin;
        Z* zkeys = _zkeys + begin;
        const int size = end - begin;
        const int BUILD_GAP_CNT = compute_gap_count(size);

        node->is_two = 0;
        node->build_size = size;
        node->size = size;
        node->fixed = 0;
        node->num_inserts = node->num_insert_to_data = 0;

        int mid1_pos = (size - 1) / 3;
        int mid2_pos = (size - 1) * 2 / 3;

        RT_ASSERT(0 <= mid1_pos);
        RT_ASSERT(mid1_pos < mid2_pos);
        RT_ASSERT(mid2_pos < size - 1);

        const long double mid1_key =
            (static_cast<long double>(zkeys[mid1_pos]) + static_cast<long double>(zkeys[mid1_pos + 1])) / 2;
        const long double mid2_key =
            (static_cast<long double>(zkeys[mid2_pos]) + static_cast<long double>(zkeys[mid2_pos + 1])) / 2;

        node->num_items = size * static_cast<int>(BUILD_GAP_CNT + 1);
        const double mid1_target = mid1_pos * static_cast<int>(BUILD_GAP_CNT + 1) + static_cast<int>(BUILD_GAP_CNT + 1) / 2;
        const double mid2_target = mid2_pos * static_cast<int>(BUILD_GAP_CNT + 1) + static_cast<int>(BUILD_GAP_CNT + 1) / 2;

        node->model.a = (mid2_target - mid1_target) / (mid2_key - mid1_key);
        node->model.b = mid1_target - node->model.a * mid1_key;
        RT_ASSERT(isfinite(node->model.a));
        RT_ASSERT(isfinite(node->model.b));

        const int lr_remains = static_cast<int>(size * BUILD_LR_REMAIN);
        node->model.b += lr_remains;
        node->num_items += lr_remains * 2;

        if (size > 1e6) {
          node->fixed = 1;
        }

        node->items = new_items(node->num_items);
        const int bitmap_size = BITMAP_SIZE(node->num_items);
        node->none_bitmap = new_bitmap(bitmap_size);
        node->child_bitmap = new_bitmap(bitmap_size);
        memset(node->none_bitmap, 0xff, sizeof(bitmap_t) * bitmap_size);
        memset(node->child_bitmap, 0, sizeof(bitmap_t) * bitmap_size);

        for (int item_i = PREDICT_POS(node, zkeys[0]), offset = 0; offset < size; ) {
          int next = offset + 1, next_i = -1;
          while (next < size) {
            next_i = PREDICT_POS(node, zkeys[next]);
            if (next_i == item_i) {
              next ++;
            } else {
              break;
            }
          }
          if (next == offset + 1) {
            BITMAP_CLEAR(node->none_bitmap, item_i);
            node->items[item_i].comp.data = keys[offset];
          } else {
            // ASSERT(next - offset <= (size+2) / 3);
            BITMAP_CLEAR(node->none_bitmap, item_i);
            BITMAP_SET(node->child_bitmap, item_i);
            node->items[item_i].comp.child = new_nodes(1);
            s.push((Segment){begin + offset, begin + next, level + 1, node->items[item_i].comp.child});
          }
          if (next >= size) {
            break;
          } else {
            item_i = next_i;
            offset = next;
          }
        }
      }
    }

    delete[] _zkeys;

    return ret;
  }
  /// bulk build, _keys must be sorted in asc order.
  /// FMCD method.
  Node* build_tree_bulk_fmcd(Point<T>* _keys, int _size, bool sorted)
  {
    Z* _zkeys = new Z[_size];
    for (int i = 0; i < _size; ++i)
      _zkeys[i] = encode(_keys[i]);
    if (!sorted)
    {
      std::pair<Z, Point<T>>* pairs = new std::pair<Z, Point<T>>[_size];
      for (int i = 0; i < _size; ++i )
        pairs[i] = std::make_pair(_zkeys[i], _keys[i]);
      std::sort(pairs, pairs + _size, [](std::pair<Z, Point<T>> &left, std::pair<Z, Point<T>> &right) {
        return left.first < right.first;
      });
      for (int i = 0; i < _size; ++i )
      {
        _zkeys[ i ] = pairs[i].first;
        _keys[ i ] = pairs[i].second;
      }
      delete[] pairs;
    }

    RT_ASSERT(_size > 1);

    typedef struct {
      int begin;
      int end;
      int level; // top level = 1
      Node* node;
    } Segment;
    std::stack<Segment> s;

    Node* ret = new_nodes(1);
    s.push((Segment){0, _size, 1, ret});

    while (!s.empty()) {
      const int begin = s.top().begin;
      const int end = s.top().end;
      const int level = s.top().level;
      Node* node = s.top().node;
      s.pop();

      RT_ASSERT(end - begin >= 2);
      if (end - begin == 2) {
        Node* _ = build_tree_two(_keys[begin], _zkeys[begin], _keys[begin+1], _zkeys[begin+1]);
        memcpy(node, _, sizeof(Node));
        delete_nodes(_, 1);
      } else {
        Point<T>* keys = _keys + begin;
        Z* zkeys = _zkeys + begin;
        const int size = end - begin;
        const int BUILD_GAP_CNT = compute_gap_count(size);

        node->is_two = 0;
        node->build_size = size;
        node->size = size;
        node->fixed = 0;
        node->num_inserts = node->num_insert_to_data = 0;

        // FMCD method
        // Here the implementation is a little different with Algorithm 1 in our paper.
        // In Algorithm 1, U_T should be (keys[size-1-D] - keys[D]) / (L - 2).
        // But according to the derivation described in our paper, M.A should be less than 1 / U_T.
        // So we added a small number (1e-6) to U_T.
        // In fact, it has only a negligible impact of the performance.
        {
          const int L = size * static_cast<int>(BUILD_GAP_CNT + 1);
          int i = 0;
          int D = 1;
          RT_ASSERT(D <= size-1-D);
          double Ut = (static_cast<long double>(zkeys[size - 1 - D]) - static_cast<long double>(zkeys[D])) /
                (static_cast<double>(L - 2)) + 1e-6;
          while (i < size - 1 - D) {
            while (i + D < size && zkeys[i + D] - zkeys[i] >= Ut) {
              i ++;
            }
            if (i + D >= size) {
              break;
            }
            D = D + 1;
            if (D * 3 > size) break;
            RT_ASSERT(D <= size-1-D);
            Ut = (static_cast<long double>(zkeys[size - 1 - D]) - static_cast<long double>(zkeys[D])) /
               (static_cast<double>(L - 2)) + 1e-6;
          }
          if (D * 3 <= size) {
            stats.fmcd_success_times ++;

            node->model.a = 1.0 / Ut;
            node->model.b = (L - node->model.a * (static_cast<long double>(zkeys[size - 1 - D]) +
                                            static_cast<long double>(zkeys[D]))) / 2;
            RT_ASSERT(isfinite(node->model.a));
            RT_ASSERT(isfinite(node->model.b));
            node->num_items = L;
          } else {
            stats.fmcd_broken_times ++;

            int mid1_pos = (size - 1) / 3;
            int mid2_pos = (size - 1) * 2 / 3;

            RT_ASSERT(0 <= mid1_pos);
            RT_ASSERT(mid1_pos < mid2_pos);
            RT_ASSERT(mid2_pos < size - 1);

            const long double mid1_key = (static_cast<long double>(zkeys[mid1_pos]) +
                                    static_cast<long double>(zkeys[mid1_pos + 1])) / 2;
            const long double mid2_key = (static_cast<long double>(zkeys[mid2_pos]) +
                                    static_cast<long double>(zkeys[mid2_pos + 1])) / 2;

            node->num_items = size * static_cast<int>(BUILD_GAP_CNT + 1);
            const double mid1_target = mid1_pos * static_cast<int>(BUILD_GAP_CNT + 1) + static_cast<int>(BUILD_GAP_CNT + 1) / 2;
            const double mid2_target = mid2_pos * static_cast<int>(BUILD_GAP_CNT + 1) + static_cast<int>(BUILD_GAP_CNT + 1) / 2;

            node->model.a = (mid2_target - mid1_target) / (mid2_key - mid1_key);
            node->model.b = mid1_target - node->model.a * mid1_key;
            RT_ASSERT(isfinite(node->model.a));
            RT_ASSERT(isfinite(node->model.b));
          }
        }
        RT_ASSERT(node->model.a >= 0);
        const int lr_remains = static_cast<int>(size * BUILD_LR_REMAIN);
        node->model.b += lr_remains;
        node->num_items += lr_remains * 2;

        if (size > 1e6) {
          node->fixed = 1;
        }

        node->items = new_items(node->num_items);
        const int bitmap_size = BITMAP_SIZE(node->num_items);
        node->none_bitmap = new_bitmap(bitmap_size);
        node->child_bitmap = new_bitmap(bitmap_size);
        memset(node->none_bitmap, 0xff, sizeof(bitmap_t) * bitmap_size);
        memset(node->child_bitmap, 0, sizeof(bitmap_t) * bitmap_size);

        for (int item_i = PREDICT_POS(node, zkeys[0]), offset = 0; offset < size; ) {
          int next = offset + 1, next_i = -1;
          while (next < size) {
            next_i = PREDICT_POS(node, zkeys[next]);
            if (next_i == item_i) {
              next ++;
            } else {
              break;
            }
          }
          if (next == offset + 1) {
            BITMAP_CLEAR(node->none_bitmap, item_i);
            node->items[item_i].comp.data = keys[offset];
          } else {
            // ASSERT(next - offset <= (size+2) / 3);
            BITMAP_CLEAR(node->none_bitmap, item_i);
            BITMAP_SET(node->child_bitmap, item_i);
            node->items[item_i].comp.child = new_nodes(1);
            s.push((Segment){begin + offset, begin + next, level + 1, node->items[item_i].comp.child});
          }
          if (next >= size) {
            break;
          } else {
            item_i = next_i;
            offset = next;
          }
        }
      }
    }

    delete[] _zkeys;

    return ret;
  }

  void destory_pending()
  {
    while (!pending_two.empty()) {
      Node* node = pending_two.top(); pending_two.pop();

      delete_items(node->items, node->num_items);
      const int bitmap_size = BITMAP_SIZE(node->num_items);
      delete_bitmap(node->none_bitmap, bitmap_size);
      delete_bitmap(node->child_bitmap, bitmap_size);
      delete_nodes(node, 1);
    }
  }

  void destroy_tree(Node* root)
  {
    std::stack<Node*> s;
    s.push(root);
    while (!s.empty()) {
      Node* node = s.top(); s.pop();

      for (int i = 0; i < node->num_items; i ++) {
        if (BITMAP_GET(node->child_bitmap, i) == 1) {
          s.push(node->items[i].comp.child);
        }
      }

      if (node->is_two) {
        RT_ASSERT(node->build_size == 2);
        RT_ASSERT(node->num_items == 8);
        node->size = 2;
        node->num_inserts = node->num_insert_to_data = 0;
        node->none_bitmap[0] = 0xff;
        node->child_bitmap[0] = 0;
        pending_two.push(node);
      } else {
        delete_items(node->items, node->num_items);
        const int bitmap_size = BITMAP_SIZE(node->num_items);
        delete_bitmap(node->none_bitmap, bitmap_size);
        delete_bitmap(node->child_bitmap, bitmap_size);
        delete_nodes(node, 1);
      }
    }
  }

  void scan_and_destory_tree(Node* _root, Point<T>* keys, bool destory = true)
  {
    typedef std::pair<int, Node*> Segment; // <begin, Node*>
    std::stack<Segment> s;

    s.push(Segment(0, _root));
    while (!s.empty()) {
      int begin = s.top().first;
      Node* node = s.top().second;
      #ifdef DEBUG
      const int SHOULD_END_POS = begin + node->size;
      #endif
      s.pop();

      for (int i = 0; i < node->num_items; i ++) {
        if (BITMAP_GET(node->none_bitmap, i) == 0) {
          if (BITMAP_GET(node->child_bitmap, i) == 0) {
            keys[begin] = node->items[i].comp.data;
            begin ++;
          } else {
            s.push(Segment(begin, node->items[i].comp.child));
            begin += node->items[i].comp.child->size;
          }
        }
      }
      RT_ASSERT(SHOULD_END_POS == begin);

      if (destory) {
        if (node->is_two) {
          RT_ASSERT(node->build_size == 2);
          RT_ASSERT(node->num_items == 8);
          node->size = 2;
          node->num_inserts = node->num_insert_to_data = 0;
          node->none_bitmap[0] = 0xff;
          node->child_bitmap[0] = 0;
          pending_two.push(node);
        } else {
          delete_items(node->items, node->num_items);
          const int bitmap_size = BITMAP_SIZE(node->num_items);
          delete_bitmap(node->none_bitmap, bitmap_size);
          delete_bitmap(node->child_bitmap, bitmap_size);
          delete_nodes(node, 1);
        }
      }
    }
  }

  Node* insert_tree(Node* _node, const Point<T>& key)
  {
    const Z zkey = encode(key);
    constexpr int MAX_DEPTH = 128;
    Node* path[MAX_DEPTH];
    int path_size = 0;
    int insert_to_data = 0;

    for (Node* node = _node; ; ) {
      RT_ASSERT(path_size < MAX_DEPTH);
      path[path_size ++] = node;

      node->size ++;
      node->num_inserts ++;
      int pos = PREDICT_POS(node, zkey);
      if (BITMAP_GET(node->none_bitmap, pos) == 1) {
        BITMAP_CLEAR(node->none_bitmap, pos);
        node->items[pos].comp.data = key;
        break;
      } else if (BITMAP_GET(node->child_bitmap, pos) == 0) {
        BITMAP_SET(node->child_bitmap, pos);
        node->items[pos].comp.child = build_tree_two(key, zkey, node->items[pos].comp.data, encode(node->items[pos].comp.data));
        insert_to_data = 1;
        break;
      } else {
        node = node->items[pos].comp.child;
      }
    }
    for (int i = 0; i < path_size; i ++) {
      path[i]->num_insert_to_data += insert_to_data;
    }

    for (int i = 0; i < path_size; i ++) {
      Node* node = path[i];
      const int num_inserts = node->num_inserts;
      const int num_insert_to_data = node->num_insert_to_data;
      const bool need_rebuild = node->fixed == 0 && node->size >= node->build_size * 4 && node->size >= 64 && num_insert_to_data * 10 >= num_inserts;

      if (need_rebuild) {
        const int ESIZE = node->size;
        Point<T>* keys = new Point<T>[ESIZE];

        scan_and_destory_tree(node, keys);
        Node* new_node = build_tree_bulk(keys, ESIZE, true);

        delete[] keys;

        path[i] = new_node;
        if (i > 0) {
          int pos = PREDICT_POS(path[i-1], zkey);
          path[i-1]->items[pos].comp.child = new_node;
        }

        break;
      }
    }

    return path[0];
  }
};

#endif // __MLIPP_Z_H__
