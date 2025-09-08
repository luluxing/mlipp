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

  bool operator!=(const Point& other) const {
    return x != other.x || y != other.y;
  }

  bool operator<(const Point& other) const {
    return x < other.x || (x == other.x && y < other.y);
  }

  bool operator>(const Point& other) const {
    return x > other.x || (x == other.x && y > other.y);
  }

  // Implement an operator to append a point to a stream
  friend std::ostream& operator<<(std::ostream& os, const Point& point) {
    os << "(" << point.x << "," << point.y << ")";
    return os;
  }
};

#define PT_VAL(point, axis) ((axis) == 0 ? (point).x : (point).y)
#define PT_EQ(p1, p2) ((p1).x == (p2).x && (p1).y == (p2).y)

#define PAIR_VAL(p, axis) (((int*)&p)[axis])

template <class T>
static inline int compare_x(const void* a, const void* b) {
  if (static_cast<const Point<T>*>(a)->x == static_cast<const Point<T>*>(b)->x)
    return 0;
  else if (static_cast<const Point<T>*>(a)->x < static_cast<const Point<T>*>(b)->x)
    return -1;
  else
    return 1;
}

template <class T>
static inline int compare_y(const void* a, const void* b) {
  if (static_cast<const Point<T>*>(a)->y == static_cast<const Point<T>*>(b)->y)
    return 0;
  else if (static_cast<const Point<T>*>(a)->y < static_cast<const Point<T>*>(b)->y)
    return -1;
  else
    return 1;
}

template <class T>
static inline int compare_pt(const void* a, const void* b) {
  int cmp_x = compare_x<T>(a, b);
  if (cmp_x == 0)
    return compare_y<T>(a, b);
  else
    return cmp_x;
}

template <class T>
static inline T get_dim(const Point<T>& key, bool x_axis) {
  return x_axis ? key.x : key.y;
}

#endif // __POINT_H_