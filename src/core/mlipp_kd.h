#ifndef __MLIPP_KD_H__
#define __MLIPP_KD_H__

#include <stdint.h>
#include <math.h>
#include <cfloat>
#include <limits>
#include <cstdio>
#include <stack>
#include <vector>
#include <cstring>
#include <sstream>

#include "lipp_base.h"
#include "lipp_utils.h"
#include "point.h"
#include "mlipp_pairingheap.h"

template<class T, bool USE_FMCD = true>
class MLIPP_KD
{

  inline int compute_gap_count(int size) {
    if (size >= 1000000) return 1;
    if (size >= 100000) return 2;
    return 5;
  }

  struct Node;
  inline int PREDICT_POS(Node* node, const Point<T>& key) const {
    double v = node->model.predict_double(PT_VAL(key, node->level % 2));
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

#ifdef BREAKDOWN
  uint64_t nodes_contain_result = 0;
  uint64_t total_checked_nodes = 0;
  double innernode_search_time = 0.0;
  double leafnode_search_time = 0.0;
  uint32_t max_depth = 0;
#endif

  MLIPP_KD(double BUILD_LR_REMAIN = 0, bool QUIET = true)
    : BUILD_LR_REMAIN(BUILD_LR_REMAIN), QUIET(QUIET) {
    /*{
      std::vector<Node*> nodes;
      for (int _ = 0; _ < 1e7; _ ++) {
        Node* node = build_tree_two(V(0, 0), V(1, 1), 0);
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
  ~MLIPP_KD() {
    destroy_tree(root);
    root = NULL;
    destory_pending();
  }

  void insert(const Point<T>& key) {
    root = insert_tree(root, key);
  }
  bool exists(const Point<T>& key) const {
    Node* node = root;
    while (true) {
      int pos = PREDICT_POS(node, key);
      if (BITMAP_GET(node->none_bitmap, pos) == 1) {
        return false;
      } else if (BITMAP_GET(node->child_bitmap, pos) == 0) {
        return PT_EQ(node->items[pos].comp.data, key);
      } else {
        node = node->items[pos].comp.child;
      }
    }
  }
  void rangeQuery(const Point<T>& min_key,
          const Point<T>& max_key,
          std::vector<Point<T>>& result) {
    rangeQueryInternal2(root, min_key, max_key, result);
    /*rangeQueryInternal(root, min_key, max_key,
      std::make_pair(true, true), 
      std::make_pair(true, true),
      result);*/
    return;
  }
  // Find the keys which are in range [lower, upper], returns the number of found keys.
  int range_query(const Point<T>& lower, const Point<T>& upper, Point<T>* results) const {
    return range_core_x<false, false, false, false>(results, 0, root, lower, upper);
  }
  int knn_query(const Point<T>& key, int k, Point<T>* results, double* distances) {
    return knn_core(root, key, k, results, distances);
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
      root = build_tree_two(vs[0], vs[1], 0);
      return;
    }

    RT_ASSERT(num_keys > 2);
    /*for (int i = 1; i < num_keys; i ++) {
      RT_ASSERT(vs[i].first > vs[i-1].first);
    }*/

    destroy_tree(root);
    root = build_tree_bulk(vs, num_keys, 0);
  }

  /*void show() const {
    printf("============= SHOW MLIPP_KD ================\n");

    std::stack<Node*> s;
    s.push(root);
    while (!s.empty()) {
      Node* node = s.top(); s.pop();

      printf("Node(%p, a = %lf, b = %llf, num_items = %d)", node, node->model.a, node->model.b, node->num_items);
      printf("[");
      int first = 1;
      for (int i = 0; i < node->num_items; i ++) {
        if (!first) {
          printf(", ");
        }
        first = 0;
        if (BITMAP_GET(node->none_bitmap, i) == 1) {
          printf("None");
        } else if (BITMAP_GET(node->child_bitmap, i) == 0) {
          std::stringstream s;
          s << node->items[i].comp.data.key;
          printf("Key(%s)", s.str().c_str());
        } else {
          printf("Child(%p)", node->items[i].comp.child);
          s.push(node->items[i].comp.child);
        }
      }
      printf("]\n");
    }
  }*/
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
    int level; // level of the node in the tree (rootlevel = 0)
    int is_two; // is special node for only two keys
    int build_size; // tree size (include sub nodes) when node created
    int size; // current tree size (include sub nodes)
    int fixed; // fixed node will not trigger rebuild
    int num_inserts, num_insert_to_data;
    int num_items; // size of items
    LinearModel<T> model;
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

  int knn_core(Node* root, const Point<T>& key, int k,
    Point<T>* results, double* distances)
  {
    double maxDist = DBL_MAX;
    int n = 0;

    int pos = PREDICT_POS(root, key);
    MPHeap heap = {NULL};
    MPNode *min_elem = pairingheap_allocate(MP_ANY, root, pos, 0, 0, 0);

    do {
      int entry_type = min_elem->entry_type;
      Node* node = (Node*)min_elem->node;
      pos = min_elem->pos;

      // printf("%lf\n", sqrt(min_elem->distance));
      /* End if distance is larger than maxDist */
      if (min_elem->distance > maxDist)
        break;

      /* Check value */
      if (BITMAP_GET(node->none_bitmap, pos) == 0) {
        if (BITMAP_GET(node->child_bitmap, pos) == 0) {
          /* It's a data point */
          Point<T> p = node->items[pos].comp.data;
          double dist = pow(p.x - key.x, 2) + pow(p.y - key.y, 2);
          // printf("Distance = %f\n", sqrt(dist));
          if (dist < maxDist) {
            /* Find its position in KNN list */
            int i = n;
            while (i > 0 && dist < distances[i-1])
              i--;
            /* Insert it into KNN list */
            int shift_count = n < k ? n - i : n - i - 1;
            memmove(&results[i + 1], &results[i], sizeof(Point<T>)*shift_count);
            memmove(&distances[i + 1], &distances[i], sizeof(double)*shift_count);
            results[i] = p;
            distances[i] = dist;
            /* Update maxDist if we replaced an element */
            if (n < k)
              n++;
            if (n == k)
              maxDist = distances[n - 1];
          }
        } else {
          /* It's a node */
          Node* child_node = node->items[pos].comp.child;
          int child_pos = PREDICT_POS(child_node, key);
          MPNode *child_elem = pairingheap_allocate(MP_ANY, child_node, child_pos, 
            min_elem->dx, min_elem->dy, min_elem->distance);
          pairingheap_add(&heap, child_elem);
        }
      }

      /* TODO: skip empty slots */

      if (entry_type & 1 && pos > 0) {
        int axis = node->level % 2;
        double dx = min_elem->dx;
        double dy = min_elem->dy;
        if (pos == node->num_items - 1 || entry_type == MP_ANY) {
          double key_pos = node->model.predict_double(PT_VAL(key, axis));
          if (axis == 0)
            dx = (key_pos - static_cast<double>(pos)) / node->model.a;
          else
            dy = (key_pos - static_cast<double>(pos)) / node->model.a;
        } else {
          if (axis == 0)
            dx += 1.0 / node->model.a;
          else
            dy += 1.0 / node->model.a;
        }
        double distance = dx*dx + dy*dy;
        MPNode *left_elem = pairingheap_allocate(MP_LEFT, node, pos - 1, 
          dx, dy, distance);
        pairingheap_add(&heap, left_elem);
      }

      if (entry_type & 2 && pos < node->num_items - 1) {
        int axis = node->level % 2;
        double dx = min_elem->dx;
        double dy = min_elem->dy;
        if (pos == 0 || entry_type == MP_ANY) {
          double key_pos = node->model.predict_double(PT_VAL(key, axis));
          if (axis == 0)
            dx = (static_cast<double>(pos + 1) - key_pos) / node->model.a;
          else
            dy = (static_cast<double>(pos + 1) - key_pos) / node->model.a;
        } else {
          if (axis == 0)
            dx += 1.0 / node->model.a;
          else
            dy += 1.0 / node->model.a;
        }
        double distance = dx*dx + dy*dy;
        MPNode *right_elem = pairingheap_allocate(MP_RIGHT, node, pos + 1, 
          dx, dy, distance);
        pairingheap_add(&heap, right_elem);
      }
      free(min_elem);
    } while (pairingheap_remove_first(&heap, &min_elem));

    /* We were working with squared distances for efficiency */
    for (int i = 0; i < n; ++i)
      distances[i] = sqrt(distances[i]);
    return n;
  }

  // SATISFY_LOWER = true means all the keys in the subtree of `node` are no less than to `lower`.
  // SATISFY_UPPER = true means all the keys in the subtree of `node` are no greater than to `upper`.
  template<bool SATISFY_LOWER_X, bool SATISFY_UPPER_X, 
       bool SATISFY_LOWER_Y, bool SATISFY_UPPER_Y>
  int range_core_x(Point<T>* results, int pos, Node* node, const Point<T>& lower, const Point<T>& upper) const
  {
    if constexpr (SATISFY_LOWER_X && SATISFY_UPPER_X) {
      int bit_pos = 0;
      const bitmap_t* none_bitmap = node->none_bitmap;
      while (bit_pos < node->num_items) {
        bitmap_t not_none = ~(*none_bitmap);
        while (not_none) {
          int latest_pos = BITMAP_NEXT_1(not_none);
          not_none ^= 1 << latest_pos;

          int i = bit_pos + latest_pos;
          if constexpr (SATISFY_LOWER_Y && SATISFY_UPPER_Y) {
            if (BITMAP_GET(node->child_bitmap, i) == 0) {
              results[pos++] = node->items[i].comp.data;
            } else {
              pos = range_core_y<true, true, true, true>(
                results, pos, node->items[i].comp.child, lower, upper);
            }
          } else if constexpr (SATISFY_LOWER_Y && !SATISFY_UPPER_Y) {
            if (BITMAP_GET(node->child_bitmap, i) == 0) {
              if (node->items[i].comp.data.y <= upper.y) {
                results[pos++] = node->items[i].comp.data;
              }
            } else {
              pos = range_core_y<true, false, true, true>(
                results, pos, node->items[i].comp.child, lower, upper);
            }
          } else if constexpr (!SATISFY_LOWER_Y && SATISFY_UPPER_Y) {
            if (BITMAP_GET(node->child_bitmap, i) == 0) {
              if (node->items[i].comp.data.y >= lower.y) {
                results[pos++] = node->items[i].comp.data;
              }
            } else {
              pos = range_core_y<false, true, true, true>(
                results, pos, node->items[i].comp.child, lower, upper);
            }
          } else { // if constexpr (!SATISFY_LOWER_Y && !SATISFY_UPPER_Y)
            if (BITMAP_GET(node->child_bitmap, i) == 0) {
              if (node->items[i].comp.data.y >= lower.y
                && node->items[i].comp.data.y <= upper.y) {
                results[pos++] = node->items[i].comp.data;
              }
            } else {
              pos = range_core_y<false, false, true, true>(
                results, pos, node->items[i].comp.child, lower, upper);
            }
          }
        }

        bit_pos += BITMAP_WIDTH;
        none_bitmap ++;
      }
      return pos;
    } else {
      int lower_pos = SATISFY_LOWER_X ? -1 : PREDICT_POS(node, lower);
      int upper_pos = SATISFY_UPPER_X ? node->num_items : PREDICT_POS(node, upper);
      if constexpr (!SATISFY_LOWER_X) {
        if (BITMAP_GET(node->none_bitmap, lower_pos) == 0) {
          if constexpr (SATISFY_LOWER_Y && SATISFY_UPPER_Y) {
            if (BITMAP_GET(node->child_bitmap, lower_pos) == 0) {
              do {
                if (node->items[lower_pos].comp.data.x < lower.x) break;
                if constexpr (!SATISFY_UPPER_X) {
                  if (node->items[lower_pos].comp.data.x > upper.x) break;
                }
                results[pos++] = node->items[lower_pos].comp.data;
              } while (false);
            } else {
              if (lower_pos < upper_pos) {
                pos = range_core_y<true, true, false, true>(
                  results, pos, node->items[lower_pos].comp.child, lower, upper);
              } else {
                pos = range_core_y<true, true, false, false>(
                  results, pos, node->items[lower_pos].comp.child, lower, upper);
              }
            }
          } else if constexpr (SATISFY_LOWER_Y && !SATISFY_UPPER_Y) {
            if (BITMAP_GET(node->child_bitmap, lower_pos) == 0) {
              do {
                if (node->items[lower_pos].comp.data.x < lower.x) break;
                if constexpr (!SATISFY_UPPER_X) {
                  if (node->items[lower_pos].comp.data.x > upper.x) break;
                }
                if (node->items[lower_pos].comp.data.y <= upper.y) {
                  results[pos++] = node->items[lower_pos].comp.data;
                }
              } while (false);
            } else {
              if (lower_pos < upper_pos) {
                pos = range_core_y<true, false, false, true>(
                  results, pos, node->items[lower_pos].comp.child, lower, upper);
              } else {
                pos = range_core_y<true, false, false, false>(
                  results, pos, node->items[lower_pos].comp.child, lower, upper);
              }
            }
          } else if constexpr (!SATISFY_LOWER_Y && SATISFY_UPPER_Y) {
            if (BITMAP_GET(node->child_bitmap, lower_pos) == 0) {
              do {
                if (node->items[lower_pos].comp.data.x < lower.x) break;
                if constexpr (!SATISFY_UPPER_X) {
                  if (node->items[lower_pos].comp.data.x > upper.x) break;
                }
                if (node->items[lower_pos].comp.data.y >= lower.y) {
                  results[pos++] = node->items[lower_pos].comp.data;
                }
              } while (false);
            } else {
              if (lower_pos < upper_pos) {
                pos = range_core_y<false, true, false, true>(
                  results, pos, node->items[lower_pos].comp.child, lower, upper);
              } else {
                pos = range_core_y<false, true, false, false>(
                  results, pos, node->items[lower_pos].comp.child, lower, upper);
              }
            }
          } else { // if constexpr (!SATISFY_LOWER_Y && !SATISFY_UPPER_Y)
            if (BITMAP_GET(node->child_bitmap, lower_pos) == 0) {
              do {
                if (node->items[lower_pos].comp.data.x < lower.x) break;
                if constexpr (!SATISFY_UPPER_X) {
                  if (node->items[lower_pos].comp.data.x > upper.x) break;
                }
                if (node->items[lower_pos].comp.data.y >= lower.y
                  && node->items[lower_pos].comp.data.y <= upper.y) {
                  results[pos++] = node->items[lower_pos].comp.data;
                }
              } while (false);
            } else {
              if (lower_pos < upper_pos) {
                pos = range_core_y<false, false, false, true>(
                  results, pos, node->items[lower_pos].comp.child, lower, upper);
              } else {
                pos = range_core_y<false, false, false, false>(
                  results, pos, node->items[lower_pos].comp.child, lower, upper);
              }
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

            if constexpr (SATISFY_LOWER_Y && SATISFY_UPPER_Y) {
              if (BITMAP_GET(node->child_bitmap, i) == 0) {
                results[pos++] = node->items[i].comp.data;
              } else {
                pos = range_core_y<true, true, true, true>(
                  results, pos, node->items[i].comp.child, lower, upper);
              }
            } else if constexpr (SATISFY_LOWER_Y && !SATISFY_UPPER_Y) {
              if (BITMAP_GET(node->child_bitmap, i) == 0) {
                if (node->items[i].comp.data.y <= upper.y) {
                  results[pos++] = node->items[i].comp.data;
                }
              } else {
                pos = range_core_y<true, false, true, true>(
                  results, pos, node->items[i].comp.child, lower, upper);
              }
            } else if constexpr (!SATISFY_LOWER_Y && SATISFY_UPPER_Y) {
              if (BITMAP_GET(node->child_bitmap, i) == 0) {
                if (node->items[i].comp.data.y >= lower.y) {
                  results[pos++] = node->items[i].comp.data;
                }
              } else {
                pos = range_core_y<false, true, true, true>(
                  results, pos, node->items[i].comp.child, lower, upper);
              }
            } else { // if constexpr (!SATISFY_LOWER_Y && !SATISFY_UPPER_Y)
              if (BITMAP_GET(node->child_bitmap, i) == 0) {
                if (node->items[i].comp.data.y >= lower.y
                  && node->items[i].comp.data.y <= upper.y) {
                  results[pos++] = node->items[i].comp.data;
                }
              } else {
                pos = range_core_y<false, false, true, true>(
                  results, pos, node->items[i].comp.child, lower, upper);
              }
            }
          }

          bit_pos += BITMAP_WIDTH;
          none_bitmap ++;
        }
      }
      if constexpr (!SATISFY_UPPER_X) {
        if (lower_pos < upper_pos) {
          if (BITMAP_GET(node->none_bitmap, upper_pos) == 0) {
            if constexpr (SATISFY_LOWER_Y && SATISFY_UPPER_Y) {
              if (BITMAP_GET(node->child_bitmap, upper_pos) == 0) {
                if (node->items[upper_pos].comp.data.x <= upper.x) {
                  results[pos++] = node->items[upper_pos].comp.data;
                }
              } else {
                pos = range_core_y<true, true, true, false>(
                  results, pos, node->items[upper_pos].comp.child, lower, upper);
              }
            } else if constexpr (SATISFY_LOWER_Y && !SATISFY_UPPER_Y) {
              if (BITMAP_GET(node->child_bitmap, upper_pos) == 0) {
                if (node->items[upper_pos].comp.data.x <= upper.x
                  && node->items[upper_pos].comp.data.y <= upper.y) {
                  results[pos++] = node->items[upper_pos].comp.data;
                }
              } else {
                pos = range_core_y<true, false, true, false>(
                  results, pos, node->items[upper_pos].comp.child, lower, upper);
              }
            } else if constexpr (!SATISFY_LOWER_Y && SATISFY_UPPER_Y) {
              if (BITMAP_GET(node->child_bitmap, upper_pos) == 0) {
                if (node->items[upper_pos].comp.data.x <= upper.x
                  && node->items[upper_pos].comp.data.y >= lower.y) {
                  results[pos++] = node->items[upper_pos].comp.data;
                }
              } else {
                pos = range_core_y<false, true, true, false>(
                  results, pos, node->items[upper_pos].comp.child, lower, upper);
              }
            } else { // if constexpr (!SATISFY_LOWER_Y && !SATISFY_UPPER_Y)
              if (BITMAP_GET(node->child_bitmap, upper_pos) == 0) {
                if (node->items[upper_pos].comp.data.x <= upper.x
                  && node->items[upper_pos].comp.data.y >= lower.y
                  && node->items[upper_pos].comp.data.y <= upper.y) {
                  results[pos++] = node->items[upper_pos].comp.data;
                }
              } else {
                pos = range_core_y<false, false, true, false>(
                  results, pos, node->items[upper_pos].comp.child, lower, upper);
              }
            }
          }
        }
      }
      return pos;
    }
  }

  // SATISFY_LOWER = true means all the keys in the subtree of `node` are no less than to `lower`.
  // SATISFY_UPPER = true means all the keys in the subtree of `node` are no greater than to `upper`.
  template<bool SATISFY_LOWER_Y, bool SATISFY_UPPER_Y, 
       bool SATISFY_LOWER_X, bool SATISFY_UPPER_X>
  int range_core_y(Point<T>* results, int pos, Node* node, const Point<T>& lower, const Point<T>& upper) const
  {
    if constexpr (SATISFY_LOWER_Y && SATISFY_UPPER_Y) {
      int bit_pos = 0;
      const bitmap_t* none_bitmap = node->none_bitmap;
      while (bit_pos < node->num_items) {
        bitmap_t not_none = ~(*none_bitmap);
        while (not_none) {
          int latest_pos = BITMAP_NEXT_1(not_none);
          not_none ^= 1 << latest_pos;

          int i = bit_pos + latest_pos;
          if constexpr (SATISFY_LOWER_X && SATISFY_UPPER_X) {
            if (BITMAP_GET(node->child_bitmap, i) == 0) {
              results[pos++] = node->items[i].comp.data;
            } else {
              pos = range_core_x<true, true, true, true>(
                results, pos, node->items[i].comp.child, lower, upper);
            }
          } else if constexpr (SATISFY_LOWER_X && !SATISFY_UPPER_X) {
            if (BITMAP_GET(node->child_bitmap, i) == 0) {
              if (node->items[i].comp.data.x <= upper.x) {
                results[pos++] = node->items[i].comp.data;
              }
            } else {
              pos = range_core_x<true, false, true, true>(
                results, pos, node->items[i].comp.child, lower, upper);
            }
          } else if constexpr (!SATISFY_LOWER_X && SATISFY_UPPER_X) {
            if (BITMAP_GET(node->child_bitmap, i) == 0) {
              if (node->items[i].comp.data.x >= lower.x) {
                results[pos++] = node->items[i].comp.data;
              }
            } else {
              pos = range_core_x<false, true, true, true>(
                results, pos, node->items[i].comp.child, lower, upper);
            }
          } else { // if constexpr (!SATISFY_LOWER_X && !SATISFY_UPPER_X)
            if (BITMAP_GET(node->child_bitmap, i) == 0) {
              if (node->items[i].comp.data.x >= lower.x
                && node->items[i].comp.data.x <= upper.x) {
                results[pos++] = node->items[i].comp.data;
              }
            } else {
              pos = range_core_x<false, false, true, true>(
                results, pos, node->items[i].comp.child, lower, upper);
            }
          }
        }

        bit_pos += BITMAP_WIDTH;
        none_bitmap ++;
      }
      return pos;
    } else {
      int lower_pos = SATISFY_LOWER_Y ? -1 : PREDICT_POS(node, lower);
      int upper_pos = SATISFY_UPPER_Y ? node->num_items : PREDICT_POS(node, upper);
      if constexpr (!SATISFY_LOWER_Y) {
        if (BITMAP_GET(node->none_bitmap, lower_pos) == 0) {
          if constexpr (SATISFY_LOWER_X && SATISFY_UPPER_X) {
            if (BITMAP_GET(node->child_bitmap, lower_pos) == 0) {
              do {
                if (node->items[lower_pos].comp.data.y < lower.y) break;
                if constexpr (!SATISFY_UPPER_Y) {
                  if (node->items[lower_pos].comp.data.y > upper.y) break;
                }
                results[pos++] = node->items[lower_pos].comp.data;
              } while (false);
            } else {
              if (lower_pos < upper_pos) {
                pos = range_core_x<true, true, false, true>(
                  results, pos, node->items[lower_pos].comp.child, lower, upper);
              } else {
                pos = range_core_x<true, true, false, false>(
                  results, pos, node->items[lower_pos].comp.child, lower, upper);
              }
            }
          } else if constexpr (SATISFY_LOWER_X && !SATISFY_UPPER_X) {
            if (BITMAP_GET(node->child_bitmap, lower_pos) == 0) {
              do {
                if (node->items[lower_pos].comp.data.y < lower.y) break;
                if constexpr (!SATISFY_UPPER_Y) {
                  if (node->items[lower_pos].comp.data.y > upper.y) break;
                }
                if (node->items[lower_pos].comp.data.x <= upper.x) {
                  results[pos++] = node->items[lower_pos].comp.data;
                }
              } while (false);
            } else {
              if (lower_pos < upper_pos) {
                pos = range_core_x<true, false, false, true>(
                  results, pos, node->items[lower_pos].comp.child, lower, upper);
              } else {
                pos = range_core_x<true, false, false, false>(
                  results, pos, node->items[lower_pos].comp.child, lower, upper);
              }
            }
          } else if constexpr (!SATISFY_LOWER_X && SATISFY_UPPER_X) {
            if (BITMAP_GET(node->child_bitmap, lower_pos) == 0) {
              do {
                if (node->items[lower_pos].comp.data.y < lower.y) break;
                if constexpr (!SATISFY_UPPER_Y) {
                  if (node->items[lower_pos].comp.data.y > upper.y) break;
                }
                if (node->items[lower_pos].comp.data.x >= lower.x) {
                  results[pos++] = node->items[lower_pos].comp.data;
                }
              } while (false);
            } else {
              if (lower_pos < upper_pos) {
                pos = range_core_x<false, true, false, true>(
                  results, pos, node->items[lower_pos].comp.child, lower, upper);
              } else {
                pos = range_core_x<false, true, false, false>(
                  results, pos, node->items[lower_pos].comp.child, lower, upper);
              }
            }
          } else { // if constexpr (!SATISFY_LOWER_X && !SATISFY_UPPER_X)
            if (BITMAP_GET(node->child_bitmap, lower_pos) == 0) {
              do {
                if (node->items[lower_pos].comp.data.y < lower.y) break;
                if constexpr (!SATISFY_UPPER_Y) {
                  if (node->items[lower_pos].comp.data.y > upper.y) break;
                }
                if (node->items[lower_pos].comp.data.x >= lower.x
                  && node->items[lower_pos].comp.data.x <= upper.x) {
                  results[pos++] = node->items[lower_pos].comp.data;
                }
              } while (false);
            } else {
              if (lower_pos < upper_pos) {
                pos = range_core_x<false, false, false, true>(
                  results, pos, node->items[lower_pos].comp.child, lower, upper);
              } else {
                pos = range_core_x<false, false, false, false>(
                  results, pos, node->items[lower_pos].comp.child, lower, upper);
              }
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

            if constexpr (SATISFY_LOWER_X && SATISFY_UPPER_X) {
              if (BITMAP_GET(node->child_bitmap, i) == 0) {
                results[pos++] = node->items[i].comp.data;
              } else {
                pos = range_core_x<true, true, true, true>(
                  results, pos, node->items[i].comp.child, lower, upper);
              }
            } else if constexpr (SATISFY_LOWER_X && !SATISFY_UPPER_X) {
              if (BITMAP_GET(node->child_bitmap, i) == 0) {
                if (node->items[i].comp.data.x <= upper.x) {
                  results[pos++] = node->items[i].comp.data;
                }
              } else {
                pos = range_core_x<true, false, true, true>(
                  results, pos, node->items[i].comp.child, lower, upper);
              }
            } else if constexpr (!SATISFY_LOWER_X && SATISFY_UPPER_X) {
              if (BITMAP_GET(node->child_bitmap, i) == 0) {
                if (node->items[i].comp.data.x >= lower.x) {
                  results[pos++] = node->items[i].comp.data;
                }
              } else {
                pos = range_core_x<false, true, true, true>(
                  results, pos, node->items[i].comp.child, lower, upper);
              }
            } else { // if constexpr (!SATISFY_LOWER_X && !SATISFY_UPPER_X)
              if (BITMAP_GET(node->child_bitmap, i) == 0) {
                if (node->items[i].comp.data.x >= lower.x
                  && node->items[i].comp.data.x <= upper.x) {
                  results[pos++] = node->items[i].comp.data;
                }
              } else {
                pos = range_core_x<false, false, true, true>(
                  results, pos, node->items[i].comp.child, lower, upper);
              }
            }
          }

          bit_pos += BITMAP_WIDTH;
          none_bitmap ++;
        }
      }
      if constexpr (!SATISFY_UPPER_Y) {
        if (lower_pos < upper_pos) {
          if (BITMAP_GET(node->none_bitmap, upper_pos) == 0) {
            if constexpr (SATISFY_LOWER_X && SATISFY_UPPER_X) {
              if (BITMAP_GET(node->child_bitmap, upper_pos) == 0) {
                if (node->items[upper_pos].comp.data.y <= upper.y) {
                  results[pos++] = node->items[upper_pos].comp.data;
                }
              } else {
                pos = range_core_x<true, true, true, false>(
                  results, pos, node->items[upper_pos].comp.child, lower, upper);
              }
            } else if constexpr (SATISFY_LOWER_X && !SATISFY_UPPER_X) {
              if (BITMAP_GET(node->child_bitmap, upper_pos) == 0) {
                if (node->items[upper_pos].comp.data.y <= upper.y
                  && node->items[upper_pos].comp.data.x <= upper.x) {
                  results[pos++] = node->items[upper_pos].comp.data;
                }
              } else {
                pos = range_core_x<true, false, true, false>(
                  results, pos, node->items[upper_pos].comp.child, lower, upper);
              }
            } else if constexpr (!SATISFY_LOWER_X && SATISFY_UPPER_X) {
              if (BITMAP_GET(node->child_bitmap, upper_pos) == 0) {
                if (node->items[upper_pos].comp.data.y <= upper.y
                  && node->items[upper_pos].comp.data.x >= lower.x) {
                  results[pos++] = node->items[upper_pos].comp.data;
                }
              } else {
                pos = range_core_x<false, true, true, false>(
                  results, pos, node->items[upper_pos].comp.child, lower, upper);
              }
            } else { // if constexpr (!SATISFY_LOWER_X && !SATISFY_UPPER_X)
              if (BITMAP_GET(node->child_bitmap, upper_pos) == 0) {
                if (node->items[upper_pos].comp.data.y <= upper.y
                  && node->items[upper_pos].comp.data.x >= lower.x
                  && node->items[upper_pos].comp.data.x <= upper.x) {
                  results[pos++] = node->items[upper_pos].comp.data;
                }
              } else {
                pos = range_core_x<false, false, true, false>(
                  results, pos, node->items[upper_pos].comp.child, lower, upper);
              }
            }
          }
        }
      }
      return pos;
    }
  }

  void rangeQueryInternal2(Node* root,
    const Point<T>& min_key, const Point<T>& max_key,
    std::vector<Point<T>>& result) {
    typedef struct {
      Node* node;
      int recheck_case[2];
    } RangeSection;
    std::stack<RangeSection> s;
    s.push((RangeSection){root, {0, 0}});
    while (!s.empty()) {
      Node* node = s.top().node;
      int recheck_case[2] = {*s.top().recheck_case};
      s.pop();

#ifdef BREAKDOWN
      total_checked_nodes++;
#endif
      int axis = node->level % 2;
      int ayis = (node->level + 1) % 2;

      int min_pos = 0;
      int max_pos = node->num_items - 1;

      switch (4*recheck_case[axis] + recheck_case[ayis]) {
        case 0:
          // axis: recheck min and max
          // ayis: recheck min and max
          min_pos = PREDICT_POS(node, min_key);
          max_pos = PREDICT_POS(node, max_key);
          for (int i = min_pos; i < max_pos + 1; ++i) {
            if (BITMAP_GET(node->none_bitmap, i) == 0) {
              if (BITMAP_GET(node->child_bitmap, i) == 0) {
                if ((i > min_pos || PT_VAL(node->items[i].comp.data, axis) >= PT_VAL(min_key, axis))
                   && (i < max_pos || PT_VAL(node->items[i].comp.data, axis) <= PT_VAL(max_key, axis))
                   && PT_VAL(node->items[i].comp.data, ayis) >= PT_VAL(min_key, ayis)
                   && PT_VAL(node->items[i].comp.data, ayis) <= PT_VAL(max_key, ayis)) {
                    result.push_back(node->items[i].comp.data);
#ifdef BREAKDOWN
                    nodes_contain_result++;
#endif
                }
              } else {
                int new_recheck_case[2] = { 0 };
                new_recheck_case[ayis] = recheck_case[ayis];
                if (i == min_pos && i == max_pos)
                  new_recheck_case[axis] = 0;
                else if (i == min_pos)
                  new_recheck_case[axis] = 1;
                else if (i == max_pos)
                  new_recheck_case[axis] = 2;
                else
                  new_recheck_case[axis] = 3;
                s.push((RangeSection){node->items[i].comp.child, {*new_recheck_case}});
              }
            }
          }
          break;
        case 1:
          // axis: recheck min and max
          // ayis: recheck min
          min_pos = PREDICT_POS(node, min_key);
          max_pos = PREDICT_POS(node, max_key);
          for (int i = min_pos; i < max_pos + 1; ++i) {
            if (BITMAP_GET(node->none_bitmap, i) == 0) {
              if (BITMAP_GET(node->child_bitmap, i) == 0) {
                if ((i > min_pos || PT_VAL(node->items[i].comp.data, axis) >= PT_VAL(min_key, axis))
                   && (i < max_pos || PT_VAL(node->items[i].comp.data, axis) <= PT_VAL(max_key, axis))
                   && PT_VAL(node->items[i].comp.data, ayis) >= PT_VAL(min_key, ayis)) {
                    result.push_back(node->items[i].comp.data);
#ifdef BREAKDOWN
                    nodes_contain_result++;
#endif
                   }
              } else {
                int new_recheck_case[2] = { 0 };
                new_recheck_case[ayis] = recheck_case[ayis];
                if (i == min_pos && i == max_pos)
                  new_recheck_case[axis] = 0;
                else if (i == min_pos)
                  new_recheck_case[axis] = 1;
                else if (i == max_pos)
                  new_recheck_case[axis] = 2;
                else
                  new_recheck_case[axis] = 3;
                s.push((RangeSection){node->items[i].comp.child, {*new_recheck_case}});
              }
            }
          }
          break;
        case 2:
          // axis: recheck min and max
          // ayis: recheck max
          min_pos = PREDICT_POS(node, min_key);
          max_pos = PREDICT_POS(node, max_key);
          for (int i = min_pos; i < max_pos + 1; ++i) {
            if (BITMAP_GET(node->none_bitmap, i) == 0) {
              if (BITMAP_GET(node->child_bitmap, i) == 0) {
                if ((i > min_pos || PT_VAL(node->items[i].comp.data, axis) >= PT_VAL(min_key, axis))
                 && (i < max_pos || PT_VAL(node->items[i].comp.data, axis) <= PT_VAL(max_key, axis))
                 && PT_VAL(node->items[i].comp.data, ayis) <= PT_VAL(max_key, ayis)) {
                  result.push_back(node->items[i].comp.data);
#ifdef BREAKDOWN
                    nodes_contain_result++;
#endif
                 }
              } else {
                int new_recheck_case[2] = { 0 };
                new_recheck_case[ayis] = recheck_case[ayis];
                if (i == min_pos && i == max_pos)
                  new_recheck_case[axis] = 0;
                else if (i == min_pos)
                  new_recheck_case[axis] = 1;
                else if (i == max_pos)
                  new_recheck_case[axis] = 2;
                else
                  new_recheck_case[axis] = 3;
                s.push((RangeSection){node->items[i].comp.child, {*new_recheck_case}});
              }
            }
          }
          break;
        case 3:
          // axis: recheck min and max
          // ayis: no recheck
          min_pos = PREDICT_POS(node, min_key);
          max_pos = PREDICT_POS(node, max_key);
          for (int i = min_pos; i < max_pos + 1; ++i) {
            if (BITMAP_GET(node->none_bitmap, i) == 0) {
              if (BITMAP_GET(node->child_bitmap, i) == 0) {
                if ((i > min_pos || PT_VAL(node->items[i].comp.data, axis) >= PT_VAL(min_key, axis))
                 && (i < max_pos || PT_VAL(node->items[i].comp.data, axis) <= PT_VAL(max_key, axis))) {
                  result.push_back(node->items[i].comp.data);
#ifdef BREAKDOWN
                  nodes_contain_result++;
#endif
                 }
              } else {
                int new_recheck_case[2] = { 0 };
                new_recheck_case[ayis] = recheck_case[ayis];
                if (i == min_pos && i == max_pos)
                  new_recheck_case[axis] = 0;
                else if (i == min_pos)
                  new_recheck_case[axis] = 1;
                else if (i == max_pos)
                  new_recheck_case[axis] = 2;
                else
                  new_recheck_case[axis] = 3;
                s.push((RangeSection){node->items[i].comp.child, {*new_recheck_case}});
              }
            }
          }
          break;
        case 4:
          // axis: recheck min
          // ayis: recheck min and max
          min_pos = PREDICT_POS(node, min_key);
          for (int i = min_pos; i < max_pos + 1; ++i) {
            if (BITMAP_GET(node->none_bitmap, i) == 0) {
              if (BITMAP_GET(node->child_bitmap, i) == 0) {
                if ((i > min_pos || PT_VAL(node->items[i].comp.data, axis) >= PT_VAL(min_key, axis))
                 && PT_VAL(node->items[i].comp.data, ayis) >= PT_VAL(min_key, ayis)
                 && PT_VAL(node->items[i].comp.data, ayis) <= PT_VAL(max_key, ayis)) {
                  result.push_back(node->items[i].comp.data);
#ifdef BREAKDOWN
                  nodes_contain_result++;
#endif
                 }
              } else {
                int new_recheck_case[2] = { 0 };
                new_recheck_case[ayis] = recheck_case[ayis];
                if (i == min_pos)
                  new_recheck_case[axis] = 1;
                else
                  new_recheck_case[axis] = 3;
                s.push((RangeSection){node->items[i].comp.child, {*new_recheck_case}});
              }
            }
          }
          break;
        case 5:
          // axis: recheck min
          // ayis: recheck min
          min_pos = PREDICT_POS(node, min_key);
          for (int i = min_pos; i < max_pos + 1; ++i) {
            if (BITMAP_GET(node->none_bitmap, i) == 0) {
              if (BITMAP_GET(node->child_bitmap, i) == 0) {
                if ((i > min_pos || PT_VAL(node->items[i].comp.data, axis) >= PT_VAL(min_key, axis))
                 && PT_VAL(node->items[i].comp.data, ayis) >= PT_VAL(min_key, ayis)) {
                  result.push_back(node->items[i].comp.data);
#ifdef BREAKDOWN
                    nodes_contain_result++;
#endif
                 }
              } else {
                int new_recheck_case[2] = { 0 };
                new_recheck_case[ayis] = recheck_case[ayis];
                if (i == min_pos)
                  new_recheck_case[axis] = 1;
                else
                  new_recheck_case[axis] = 3;
                s.push((RangeSection){node->items[i].comp.child, {*new_recheck_case}});
              }
            }
          }
          break;
        case 6:
          // axis: recheck min
          // ayis: recheck max
          min_pos = PREDICT_POS(node, min_key);
          for (int i = min_pos; i < max_pos + 1; ++i) {
            if (BITMAP_GET(node->none_bitmap, i) == 0) {
              if (BITMAP_GET(node->child_bitmap, i) == 0) {
                if ((i > min_pos || PT_VAL(node->items[i].comp.data, axis) >= PT_VAL(min_key, axis))
                 && PT_VAL(node->items[i].comp.data, ayis) <= PT_VAL(max_key, ayis)) {
                  result.push_back(node->items[i].comp.data);
#ifdef BREAKDOWN
                  nodes_contain_result++;
#endif
                 }
              } else {
                int new_recheck_case[2] = { 0 };
                new_recheck_case[ayis] = recheck_case[ayis];
                if (i == min_pos)
                  new_recheck_case[axis] = 1;
                else
                  new_recheck_case[axis] = 3;
                s.push((RangeSection){node->items[i].comp.child, {*new_recheck_case}});
              }
            }
          }
          break;
        case 7:
          // axis: recheck min
          // ayis: no recheck
          min_pos = PREDICT_POS(node, min_key);
          for (int i = min_pos; i < max_pos + 1; ++i) {
            if (BITMAP_GET(node->none_bitmap, i) == 0) {
              if (BITMAP_GET(node->child_bitmap, i) == 0) {
                if ((i > min_pos || PT_VAL(node->items[i].comp.data, axis) >= PT_VAL(min_key, axis))) {
                  result.push_back(node->items[i].comp.data);
  #ifdef BREAKDOWN
                    nodes_contain_result++;
#endif
                }
              } else {
                int new_recheck_case[2] = { 0 };
                new_recheck_case[ayis] = recheck_case[ayis];
                if (i == min_pos)
                  new_recheck_case[axis] = 1;
                else
                  new_recheck_case[axis] = 3;
                s.push((RangeSection){node->items[i].comp.child, {*new_recheck_case}});
              }
            }
          }
          break;
        case 8:
          // axis: recheck max
          // ayis: recheck min and max
          max_pos = PREDICT_POS(node, max_key);
          for (int i = min_pos; i < max_pos + 1; ++i) {
            if (BITMAP_GET(node->none_bitmap, i) == 0) {
              if (BITMAP_GET(node->child_bitmap, i) == 0) {
                if ((i < max_pos || PT_VAL(node->items[i].comp.data, axis) <= PT_VAL(max_key, axis))
                 && PT_VAL(node->items[i].comp.data, ayis) >= PT_VAL(min_key, ayis)
                 && PT_VAL(node->items[i].comp.data, ayis) <= PT_VAL(max_key, ayis)) {
                  result.push_back(node->items[i].comp.data);
#ifdef BREAKDOWN
                    nodes_contain_result++;
#endif
                 }
              } else {
                int new_recheck_case[2] = { 0 };
                new_recheck_case[ayis] = recheck_case[ayis];
                if (i == max_pos)
                  new_recheck_case[axis] = 2;
                else
                  new_recheck_case[axis] = 3;
                s.push((RangeSection){node->items[i].comp.child, {*new_recheck_case}});
              }
            }
          }
          break;
        case 9:
          // axis: recheck max
          // ayis: recheck min
          max_pos = PREDICT_POS(node, max_key);
          for (int i = min_pos; i < max_pos + 1; ++i) {
            if (BITMAP_GET(node->none_bitmap, i) == 0) {
              if (BITMAP_GET(node->child_bitmap, i) == 0) {
                if ((i < max_pos || PT_VAL(node->items[i].comp.data, axis) <= PT_VAL(max_key, axis))
                 && PT_VAL(node->items[i].comp.data, ayis) >= PT_VAL(min_key, ayis)) {
                  result.push_back(node->items[i].comp.data);
#ifdef BREAKDOWN
                    nodes_contain_result++;
#endif
                 }
              } else {
                int new_recheck_case[2] = { 0 };
                new_recheck_case[ayis] = recheck_case[ayis];
                if (i == max_pos)
                  new_recheck_case[axis] = 2;
                else
                  new_recheck_case[axis] = 3;
                s.push((RangeSection){node->items[i].comp.child, {*new_recheck_case}});
              }
            }
          }
          break;
        case 10:
          // axis: recheck max
          // ayis: recheck max
          max_pos = PREDICT_POS(node, max_key);
          for (int i = min_pos; i < max_pos + 1; ++i) {
            if (BITMAP_GET(node->none_bitmap, i) == 0) {
              if (BITMAP_GET(node->child_bitmap, i) == 0) {
                if ((i < max_pos || PT_VAL(node->items[i].comp.data, axis) <= PT_VAL(max_key, axis))
                 && PT_VAL(node->items[i].comp.data, ayis) <= PT_VAL(max_key, ayis)) {
                  result.push_back(node->items[i].comp.data);
#ifdef BREAKDOWN
                    nodes_contain_result++;
#endif
                 }
              } else {
                int new_recheck_case[2] = { 0 };
                new_recheck_case[ayis] = recheck_case[ayis];
                if (i == max_pos)
                  new_recheck_case[axis] = 2;
                else
                  new_recheck_case[axis] = 3;
                s.push((RangeSection){node->items[i].comp.child, {*new_recheck_case}});
              }
            }
          }
          break;
        case 11:
          // axis: recheck max
          // ayis: no recheck
          max_pos = PREDICT_POS(node, max_key);
          for (int i = min_pos; i < max_pos + 1; ++i) {
            if (BITMAP_GET(node->none_bitmap, i) == 0) {
              if (BITMAP_GET(node->child_bitmap, i) == 0) {
                if ((i < max_pos || PT_VAL(node->items[i].comp.data, axis) <= PT_VAL(max_key, axis))) {
                  result.push_back(node->items[i].comp.data);
#ifdef BREAKDOWN
                    nodes_contain_result++;
#endif
                }
              } else {
                int new_recheck_case[2] = { 0 };
                new_recheck_case[ayis] = recheck_case[ayis];
                if (i == max_pos)
                  new_recheck_case[axis] = 2;
                else
                  new_recheck_case[axis] = 3;
                s.push((RangeSection){node->items[i].comp.child, {*new_recheck_case}});
              }
            }
          }
          break;
        case 12:
          // axis: no recheck
          // ayis: recheck min and max
          for (int i = min_pos; i < max_pos + 1; ++i) {
            if (BITMAP_GET(node->none_bitmap, i) == 0) {
              if (BITMAP_GET(node->child_bitmap, i) == 0) {
                if (PT_VAL(node->items[i].comp.data, ayis) >= PT_VAL(min_key, ayis)
                 && PT_VAL(node->items[i].comp.data, ayis) <= PT_VAL(max_key, ayis)) {
                  result.push_back(node->items[i].comp.data);
#ifdef BREAKDOWN
                    nodes_contain_result++;
#endif
                 }
              } else {
                s.push((RangeSection){node->items[i].comp.child, {*recheck_case}});
              }
            }
          }
          break;
        case 13:
          // axis: no recheck
          // ayis: recheck min
          for (int i = min_pos; i < max_pos + 1; ++i) {
            if (BITMAP_GET(node->none_bitmap, i) == 0) {
              if (BITMAP_GET(node->child_bitmap, i) == 0) {
                if (PT_VAL(node->items[i].comp.data, ayis) >= PT_VAL(min_key, ayis)) {
                  result.push_back(node->items[i].comp.data);
#ifdef BREAKDOWN
                    nodes_contain_result++;
#endif
                }
              } else {
                s.push((RangeSection){node->items[i].comp.child, {*recheck_case}});
              }
            }
          }
          break;
        case 14:
          // axis: no recheck
          // ayis: recheck max
          for (int i = min_pos; i < max_pos + 1; ++i) {
            if (BITMAP_GET(node->none_bitmap, i) == 0) {
              if (BITMAP_GET(node->child_bitmap, i) == 0) {
                if (PT_VAL(node->items[i].comp.data, ayis) <= PT_VAL(max_key, ayis)) {
                  result.push_back(node->items[i].comp.data);
#ifdef BREAKDOWN
                  nodes_contain_result++;
#endif
                }
              } else {
                s.push((RangeSection){node->items[i].comp.child, {*recheck_case}});
              }
            }
          }
          break;
        default:
          // axis: no recheck
          // ayis: no recheck
          for (int i = min_pos; i < max_pos + 1; ++i) {
            if (BITMAP_GET(node->none_bitmap, i) == 0) {
              if (BITMAP_GET(node->child_bitmap, i) == 0) {
                result.push_back(node->items[i].comp.data);
#ifdef BREAKDOWN
                nodes_contain_result++;
#endif
              } else {
                s.push((RangeSection){node->items[i].comp.child, {*recheck_case}});
              }
            }
          }
          break;
      }
    }
    return;
  }

  void rangeQueryInternal3(Node* root,
    const Point<T>& min_key, const Point<T>& max_key,
    std::vector<Point<T>> &result) {
    
    int min_pos = PREDICT_POS(root, min_key);
    int max_pos = PREDICT_POS(root, max_key);
    for (int i = min_pos; i < max_pos + 1; i ++) {
      if (BITMAP_GET(root->none_bitmap, i) == 0) {
        if (BITMAP_GET(root->child_bitmap, i) == 0) {
          if (root->items[i].comp.data.x >= min_key.x 
           && root->items[i].comp.data.x <= max_key.x
           && root->items[i].comp.data.y >= min_key.y
           && root->items[i].comp.data.y <= max_key.y)
            result.push_back(root->items[i].comp.data);
        } else {
          rangeQueryInternal3(root->items[i].comp.child,
                     min_key, max_key,
                     result);
        }
      }
    }
    return;
  }

  void rangeQueryInternal(Node* root,
    const Point<T>& min_key, const Point<T>& max_key,
    std::pair<bool, bool> recheck_min, 
    std::pair<bool, bool> recheck_max, 
    std::vector<Point<T>> &result) {

    if (!recheck_min.first 
     && !recheck_min.second
     && !recheck_max.first
     && !recheck_max.second)
      return fullScanInternal(root, result);
    
    int min_pos = 0;
    int max_pos = root->num_items-1;
    if (PAIR_VAL(recheck_min, root->level))
      min_pos = PREDICT_POS(root, min_key);
    if (PAIR_VAL(recheck_max, root->level))
      max_pos = PREDICT_POS(root, max_key);
    for (int i = min_pos; i < max_pos + 1; i ++) {
      if (BITMAP_GET(root->none_bitmap, i) == 0) {
        if (root->level % 2 == 0)
        {
          recheck_min.first = i == min_pos;
          recheck_max.first = i == max_pos;
        }
        else
        {
          recheck_min.second = i == min_pos;
          recheck_max.second = i == max_pos;
        }
        if (BITMAP_GET(root->child_bitmap, i) == 0) {
          if (root->items[i].comp.data.x >= min_key.x 
           && root->items[i].comp.data.x <= max_key.x
           && root->items[i].comp.data.y >= min_key.y
           && root->items[i].comp.data.y <= max_key.y)
            result.push_back(root->items[i].comp.data);
        } else {
          rangeQueryInternal(root->items[i].comp.child,
                     min_key, max_key, 
                     recheck_min, recheck_max, 
                     result);
        }
      }
    }
    return;
  }

  void fullScanInternal(Node* root, std::vector<Point<T>> &result) {
    for (int i = 0; i < root->num_items; i ++) {
      if (BITMAP_GET(root->none_bitmap, i) == 0) {
        if (BITMAP_GET(root->child_bitmap, i) == 0) {
          result.push_back(root->items[i].comp.data);
        } else {
          fullScanInternal(root->items[i].comp.child, result);
        }
      }
    }
    return;
  }


  /// build an empty tree
  Node* build_tree_none()
  {
    Node* node = new_nodes(1);
    node->level = 0;
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
  Node* build_tree_two(Point<T> key1, Point<T> key2, int level)
  {
    int axis = level % 2;

    if (PT_VAL(key1, axis) == PT_VAL(key2, axis))
    {
      level += 1;
      axis = level % 2;
    }
    if (PT_VAL(key1, axis) > PT_VAL(key2, axis)) {
      std::swap(key1, key2);
    }
    // if (PT_VAL(key1, axis) == PT_VAL(key2, axis))
    //   printf("(%d, %d), (%d, %d), %d\n", key1.x, key1.y, key2.x, key2.y, axis);
    RT_ASSERT(PT_VAL(key1, axis) < PT_VAL(key2, axis));
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
    node->level = level;

    const long double mid1_val = PT_VAL(key1, axis);
    const long double mid2_val = PT_VAL(key2, axis);

    const double mid1_target = node->num_items / 3;
    const double mid2_target = node->num_items * 2 / 3;

    node->model.a = (mid2_target - mid1_target) / (mid2_val - mid1_val);
    node->model.b = mid1_target - node->model.a * mid1_val;
    RT_ASSERT(isfinite(node->model.a));
    RT_ASSERT(isfinite(node->model.b));

    { // insert key1
      int pos = PREDICT_POS(node, key1);
      RT_ASSERT(BITMAP_GET(node->none_bitmap, pos) == 1);
      BITMAP_CLEAR(node->none_bitmap, pos);
      node->items[pos].comp.data = key1;
    }
    { // insert key2
      int pos = PREDICT_POS(node, key2);
      RT_ASSERT(BITMAP_GET(node->none_bitmap, pos) == 1);
      BITMAP_CLEAR(node->none_bitmap, pos);
      node->items[pos].comp.data = key2;
    }

    return node;
  }
  /// bulk build, _keys must be sorted in asc order.
  Node* build_tree_bulk(Point<T>* _keys, int _size, int _level)
  {
    if (USE_FMCD) {
      return build_tree_bulk_fmcd(_keys, _size, _level);
    } else {
      return build_tree_bulk_fast(_keys, _size, _level);
    }
  }
  /// bulk build, _keys must be sorted in asc order.
  /// split keys into three parts at each node.
  Node* build_tree_bulk_fast(Point<T>* _keys, int _size, int _level)
  {
    RT_ASSERT(_size > 1);

    typedef struct {
      int begin;
      int end;
      int level;
      Node* node;
    } Segment;
    std::stack<Segment> s;

    Node* ret = new_nodes(1);
    s.push((Segment){0, _size, _level, ret});

    while (!s.empty()) {
      const int begin = s.top().begin;
      const int end = s.top().end;
      const int level = s.top().level;
      Node* node = s.top().node;
      s.pop();

      RT_ASSERT(end - begin >= 2);
      if (end - begin == 2) {
        Node* _ = build_tree_two(_keys[begin], _keys[begin+1], level);
        memcpy(node, _, sizeof(Node));
        delete_nodes(_, 1);
      } else {
        Point<T>* keys = _keys + begin;
        const int axis = level % 2;
        const int size = end - begin;
        const int BUILD_GAP_CNT = compute_gap_count(size);
        if (axis == 0)
          qsort(keys, size, sizeof(Point<T>), &compare_x<T>);
        else
          qsort(keys, size, sizeof(Point<T>), &compare_y<T>);

        node->level = level;
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
            (static_cast<long double>(PT_VAL(keys[mid1_pos], axis))
              + static_cast<long double>(PT_VAL(keys[mid1_pos + 1], axis))) / 2;
        const long double mid2_key =
            (static_cast<long double>(PT_VAL(keys[mid2_pos], axis)) 
              + static_cast<long double>(PT_VAL(keys[mid2_pos + 1], axis))) / 2;

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
        memset(node->none_bitmap, BITMAP_ONE, sizeof(bitmap_t) * bitmap_size);
        memset(node->child_bitmap, BITMAP_ZERO, sizeof(bitmap_t) * bitmap_size);

        for (int item_i = PREDICT_POS(node, keys[0]), offset = 0; offset < size; ) {
          int next = offset + 1, next_i = -1;
          while (next < size) {
            next_i = PREDICT_POS(node, keys[next]);
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

    return ret;
  }
  /// bulk build, _keys must be sorted in asc order.
  /// FMCD method.
  Node* build_tree_bulk_fmcd(Point<T>* _keys, int _size, int _level)
  {
    RT_ASSERT(_size > 1);

    typedef struct {
      int begin;
      int end;
      int level;
      Node* node;
    } Segment;
    std::stack<Segment> s;

    Node* ret = new_nodes(1);
    s.push((Segment){0, _size, _level, ret});

    while (!s.empty()) {
      const int begin = s.top().begin;
      const int end = s.top().end;
      const int level = s.top().level;
      Node* node = s.top().node;
      s.pop();

      RT_ASSERT(end - begin >= 2);
      if (end - begin == 2) {
        Node* _ = build_tree_two(_keys[begin], _keys[begin+1], level);
        memcpy(node, _, sizeof(Node));
        delete_nodes(_, 1);
      } else {
        Point<T>* keys = _keys + begin;
        const int axis = level % 2;
        const int size = end - begin;
        const int BUILD_GAP_CNT = compute_gap_count(size);
        if (axis == 0)
          qsort(keys, size, sizeof(Point<T>), &compare_x<T>);
        else
          qsort(keys, size, sizeof(Point<T>), &compare_y<T>);

        node->level = level;
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
          double Ut = (static_cast<long double>(PT_VAL(keys[size - 1 - D], axis))  
            - static_cast<long double>(PT_VAL(keys[D], axis))) / (static_cast<double>(L - 2)) + 1e-6;
          while (i < size - 1 - D) {
            while (i + D < size && PT_VAL(keys[i + D], axis) - PT_VAL(keys[i], axis) >= Ut) {
              i ++;
            }
            if (i + D >= size) {
              break;
            }
            D = D + 1;
            if (D * 3 > size) break;
            RT_ASSERT(D <= size-1-D);
            Ut = (static_cast<long double>(PT_VAL(keys[size - 1 - D], axis)) 
              - static_cast<long double>(PT_VAL(keys[D], axis))) / (static_cast<double>(L - 2)) + 1e-6;
          }
          if (D * 3 <= size) {
            stats.fmcd_success_times ++;

            node->model.a = 1.0 / Ut;
            node->model.b = (L - node->model.a * (static_cast<long double>(PT_VAL(keys[size - 1 - D], axis)) +
                                static_cast<long double>(PT_VAL(keys[D], axis)))) / 2;
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

            const long double mid1_key = (static_cast<long double>(PT_VAL(keys[mid1_pos], axis)) +
                            static_cast<long double>(PT_VAL(keys[mid1_pos + 1], axis))) / 2;
            const long double mid2_key = (static_cast<long double>(PT_VAL(keys[mid2_pos], axis)) +
                            static_cast<long double>(PT_VAL(keys[mid2_pos + 1], axis))) / 2;

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
        memset(node->none_bitmap, BITMAP_ONE, sizeof(bitmap_t) * bitmap_size);
        memset(node->child_bitmap, BITMAP_ZERO, sizeof(bitmap_t) * bitmap_size);

        for (int item_i = PREDICT_POS(node, keys[0]), offset = 0; offset < size; ) {
          int next = offset + 1, next_i = -1;
          while (next < size) {
            next_i = PREDICT_POS(node, keys[next]);
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
        node->none_bitmap[0] = BITMAP_ONE;
        node->child_bitmap[0] = BITMAP_ZERO;
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
          node->none_bitmap[0] = BITMAP_ONE;
          node->child_bitmap[0] = BITMAP_ZERO;
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
    constexpr int MAX_DEPTH = 128;
    Node* path[MAX_DEPTH];
    int path_size = 0;
    int insert_to_data = 0;

    for (Node* node = _node; ; ) {
      RT_ASSERT(path_size < MAX_DEPTH);
      path[path_size ++] = node;

      node->size ++;
      node->num_inserts ++;
      int pos = PREDICT_POS(node, key);
      if (BITMAP_GET(node->none_bitmap, pos) == 1) {
        BITMAP_CLEAR(node->none_bitmap, pos);
        node->items[pos].comp.data = key;
        break;
      } else if (BITMAP_GET(node->child_bitmap, pos) == 0) {
        if (PT_EQ(node->items[pos].comp.data, key))
          return path[0];
        BITMAP_SET(node->child_bitmap, pos);
        node->items[pos].comp.child = build_tree_two(key, node->items[pos].comp.data, node->level + 1);
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
        Node* new_node = build_tree_bulk(keys, ESIZE, node->level);

        delete[] keys;

        path[i] = new_node;
        if (i > 0) {
          int pos = PREDICT_POS(path[i-1], key);
          path[i-1]->items[pos].comp.child = new_node;
        }

        break;
      }
    }

    return path[0];
  }
};

#endif // __MLIPP_KD_H__
