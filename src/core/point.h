#ifndef __POINT_H_
#define __POINT_H_

// Define a structure to represent a point in 2D space
template <class T>
struct Point {
  T x;
  T y;

  bool operator==(const Point& other) const {
    return x == other.x && y == other.y;
  }

  static inline int compare_x(const void* a, const void* b) {
    if (static_cast<const Point*>(a)->x == static_cast<const Point*>(b)->x)
      return 0;
    else if (static_cast<const Point*>(a)->x < static_cast<const Point*>(b)->x)
      return -1;
    else
      return 1;
  }

  static inline int compare_y(const void* a, const void* b) {
    if (static_cast<const Point*>(a)->y == static_cast<const Point*>(b)->y)
      return 0;
    else if (static_cast<const Point*>(a)->y < static_cast<const Point*>(b)->y)
      return -1;
    else
      return 1;
  }

  static inline int compare_pt(const void* a, const void* b) {
    int cmp_x = compare_x(a, b);
    if (cmp_x == 0)
      return compare_y(a, b);
    else
      return cmp_x;
  }
};

#define PT_VAL(point, axis) ((axis) == 0 ? (point).x : (point).y)
#define PT_EQ(p1, p2) ((p1).x == (p2).x && (p1).y == (p2).y)


#endif // __POINT_H_