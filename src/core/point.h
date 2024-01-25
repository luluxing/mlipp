#ifndef __POINT_H_
#define __POINT_H_

// Define a structure to represent a point in 2D space
typedef struct Point {
    int x;
    int y;
} Point;

#define PT_VAL(point, axis) (((int*)&point)[axis])
#define PT_EQ(p1, p2) ((p1).x == (p2).x && (p1).y == (p2.y))

int compare_x(const void* a, const void* b){
    return ((Point*)a)->x - ((Point*)b)->x;
}

int compare_y(const void* a, const void* b){
    return ((Point*)a)->y - ((Point*)b)->y;
}

int compare_pt(const void* a, const void* b){
    int cmp_x = compare_x(a, b);
    if (cmp_x == 0)
        return compare_y(a, b);
    else
        return cmp_x;
}

#endif // __POINT_H_