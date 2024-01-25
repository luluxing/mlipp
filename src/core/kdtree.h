#ifndef __KDTREE_H_
#define __KDTREE_H_

#include "point.h"
// #include "pairingheap.h"

// Define a structure to represent a node in a kd-tree
typedef struct KDNode {
    Point point;
    struct KDNode* left;
    struct KDNode* right;
} KDNode;

// Define a function to create a new kd-tree node
KDNode* 
new_kd_node(Point point) {
    KDNode* node = (KDNode*)malloc(sizeof(KDNode));
    node->point = point;
    node->left = NULL;
    node->right = NULL;
    return node;
}

// Define a function to create a new kd-tree node
void 
set_kd_node(KDNode* node, Point point) {
    node->point = point;
    node->left = NULL;
    node->right = NULL;
    return;
}

int
kd_bulk_internal(Point* points, int size, int depth, int *max_depth, KDNode* root, int pos) {
    if (size <= 0)
        return pos;
    else if (size == 1)
    {
        set_kd_node(&root[pos++], points[0]);
        if (depth > (*max_depth))
            *max_depth = depth;
        return pos;
    }
    else
    {
        // Choose the axis based on depth
        int axis = depth % 2;

        // Sort points based on axis
        // if (axis == 0)
        //     qsort(points, size, sizeof(Point), &compare_x);
        // else
        //     qsort(points, size, sizeof(Point), &compare_y);

        // Find middle of array
        int middle = size / 2;
        while (middle > 0 
            && PT_VAL(points[middle-1], axis) == PT_VAL(points[middle], axis))
            middle--;

        KDNode *node = &root[pos++];
        set_kd_node(node, points[middle]);
        if (middle > 0)
        {
            node->left = &root[pos];
            pos = kd_bulk_internal(points, 
                middle, depth + 1, max_depth, root, pos);
        }
        if (size - middle - 1 > 0)
        {
            node->right = &root[pos];
            pos = kd_bulk_internal(points + middle + 1, 
                size - middle - 1, depth + 1, max_depth, root, pos);
        }
        return pos;
    }
}

// Define a function to bulk-load points into a kd-tree
KDNode*
kd_bulk2(Point* points, int size, int depth, int *max_depth) {
    if (size <= 0)
        return NULL;
    KDNode* root = (KDNode*)malloc(sizeof(KDNode)*size);
    kd_bulk_internal(points, size, depth, max_depth, root, 0);
    return root;
}

// Define a function to bulk-load points into a kd-tree
KDNode*
kd_bulk(Point* points, int size, int depth, int *max_depth) {
    if (size <= 0)
        return NULL;
    else if (size == 1)
    {
        if (depth > (*max_depth))
            *max_depth = depth;
        return new_kd_node(points[0]);
    }
    else
    {
        // Choose the axis based on depth
        int axis = depth % 2;

        // Sort points based on axis
        if (axis == 0)
            qsort(points, size, sizeof(Point), &compare_x);
        else
            qsort(points, size, sizeof(Point), &compare_y);

        // Find middle of array
        int middle = size / 2;
        while (middle > 0 
            && PT_VAL(points[middle-1], axis) == PT_VAL(points[middle], axis))
            middle--;

        KDNode *root = new_kd_node(points[middle]);
        if (middle > 0)
            root->left = kd_bulk(points, middle, depth + 1, max_depth);
        if (size - middle - 1 > 0)
            root->right = kd_bulk(points + middle + 1, size - middle - 1, depth + 1, max_depth);
        return root;
    }
}

// Define a function to insert a new point into a kd-tree
KDNode* 
kd_insert(KDNode* root, Point point, int depth) {
    if (root == NULL) {
        return new_kd_node(point);
    }

    // Do not insert duplicate tuples
    if (root->point.x == point.x && root->point.y == point.y)
        return root;

    // Choose the axis based on depth
    int axis = depth % 2;

    // Compare the point to the current node and recursively insert it
    if (PT_VAL(point, axis) < PT_VAL(root->point, axis))
        root->left = kd_insert(root->left, point, depth+1);
    else
        root->right = kd_insert(root->right, point, depth+1);

    return root;
}

// Define a function to search for a point in a kd-tree
KDNode* 
kd_search(KDNode* root, Point point, int depth) {
    if (root == NULL || (root->point.x == point.x && root->point.y == point.y))
        return root;

    // Choose the axis based on depth
    int axis = depth % 2;

    // Compare the point to the current node and recursively insert it
    if (PT_VAL(point, axis) < PT_VAL(root->point, axis))
        return kd_search(root->left, point, depth+1);
    else
        return kd_search(root->right, point, depth+1);
}

// Define a function to search for a point in a kd-tree
inline bool
kd_exists(KDNode* root, Point point, int depth) {
    return kd_search(root, point, depth) != NULL;
}

Point *
kd_range(KDNode* root, Point min_point, Point max_point, 
    int depth, int max_depth, int *count)
{
    int k = 0, size = 256, state = 0;
    Point *result = (Point*)malloc(sizeof(Point) * size);
    KDNode* node = root;
    KDNode** parents = (KDNode**)malloc(sizeof(KDNode*)*max_depth);
    while (state < 4)
    {
        // Choose the axis based on depth
        int axis = depth % 2;

        switch (state)
        {
            case 0:
                // Go to left child
                if (node->left == NULL)
                    state++;
                else
                {
                    if (PT_VAL(min_point, axis) < PT_VAL(node->point, axis))
                    {
                        parents[depth++] = node;
                        node = node->left;
                    }
                    else
                        state++;
                }
                break;
            case 1:
                // Check point
                if (node->point.x >= min_point.x
                 && node->point.y >= min_point.y
                 && node->point.x <= max_point.x
                 && node->point.y <= max_point.y)
                {
                    if (k >= size)
                    {
                        size *= 2;
                        result = (Point *)realloc(result, sizeof(Point) * size);
                    }
                    result[k++] = node->point;
                }
                state++;
                break;
            case 2:
                // Go to right child
                if (node->right == NULL)
                    state++;
                else
                {
                    if (PT_VAL(max_point, axis) >= PT_VAL(node->point, axis))
                    {
                        parents[depth++] = node;
                        node = node->right;
                        state = 0;
                    }
                    else
                        state++;
                }
                break;
            case 3:
                // Search is over
                if (node == root)
                    state++;
                else
                {
                    // Go to parent
                    KDNode *curr = node;
                    node = parents[--depth];
                    if (node->left == curr)
                        state = 1;
                }
                break;
            default:
                // Should never happen
                break;
        }
    }

    result = (Point *)realloc(result, sizeof(Point) * k);
    *count = k;
    return result;
}

/*
 * Returns the k nearest points to 'point' ordered by distance
 *
 * Note: Requires that k <= number of points in the tree
 */
// Point *
// kd_knn(KDNode *root, Point point, int k, int depth, double *distances)
// {
//     Point *result = (Point *)malloc(sizeof(Point)*k);
//     distance = (double *)malloc(sizeof(double)*k);
//     double maxdist = -1;

//     PNode *root = pairingheap_allocate(root, depth, 0, 0, 0);
//     PNode *min_elem = NULL;

//     while (pairingheap_remove_first(root, min_elem))
//     {
//         KDNode *node = min_elem->node;
//         int depth = min_elem->depth;
//         double dx = min_elem->dx;
//         double dy = min_elem->dy;
//         double distance = min_elem->distance;

//         if (mindist != -1 && distance > mindist)
//             break;

//         int axis = depth % 2;

//         double daxis = fabs(PT_VAL(point, axis) - PT_VAL(node->point, axis));

//         if (PT_VAL(point, axis) < PT_VAL(node->point, axis))
//         {
//             if (node->left != NULL)
//             {
//                 PNode *left = pairingheap_allocate(node->left, depth + 1, dx, dy, distance);
//                 pairingheap_insert(root, left);
//             }
//             if (node->right != NULL)
//             {
//                 PNode *left = pairingheap_allocate(node->left, depth + 1, dx, dy, distance);
//                 pairingheap_insert(root, left);
//             }
//         }
//         else
//         {
//             node->left != NULL
//         }
//     }
//     if (maxdist == -1)
//         printf("Problem\n");
//     return result;
// }

void
kd_depth(KDNode* root, int depth, 
    int* max_depth, long int* sum_depth, int *sum_nodes)
{
    if (root == NULL)
    {
        *max_depth = depth;
        *sum_depth = 1;
        *sum_nodes = 1;
        return;
    }
    int mdl = depth, mdr = depth;
    long int sdl = 0, sdr = 0;
    int snl = 0, snr = 0;
    if (root->left != NULL)
        kd_depth(root->left, depth + 1, &mdl, &sdl, &snl);
    if (root->right != NULL)
        kd_depth(root->right, depth + 1, &mdr, &sdr, &snr);
    *max_depth = std::max(mdl, mdr);
    *sum_depth = sdl + sdr + depth;
    *sum_nodes = snl + snr + 1;
    return;
}

// Define a function to destroy a kd-tree
void
kd_destroy(KDNode* root) {
    if (root == NULL)
        return;

    // Free children before freeing node
    kd_destroy(root->left);
    kd_destroy(root->right);
    free(root);
    return;
}

// Define a function to print the points in a kd-tree
void 
kd_print(KDNode* root, int depth) {
    if (root != NULL) {
        kd_print(root->left, depth + 1);
        printf("%*s(%d, %d)\n", depth*2, "", root->point.x, root->point.y);
        kd_print(root->right, depth + 1);
    }
}

#endif // __KDTREE_H_