// Define a structure to represent a point in 2D space
typedef struct Point {
    int x;
    int y;
} Point;

int compare_x(const void* a, const void* b){
    return ((Point*)a)->x - ((Point*)b)->x;
}

int compare_y(const void* a, const void* b){
    return ((Point*)a)->y - ((Point*)b)->y;
}

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

// Define a function to bulk-load points into a kd-tree
KDNode*
kd_bulk(Point* points, int size, int depth) {
    if (size <= 0)
        return NULL;
    else if (size == 1)
        return new_kd_node(points[0]);
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
        while (middle > 0 && 
            ((axis == 0 && points[middle - 1].x == points[middle].x) 
                || (axis == 1 && points[middle - 1].y == points[middle].y)))
            middle--;

        KDNode *root = new_kd_node(points[middle]);
        if (middle > 0)
            root->left = kd_bulk(points, middle, depth + 1);
        if (size - middle - 1 > 0)
            root->right = kd_bulk(points + middle + 1, size - middle - 1, depth + 1);
        return root;
    }
}

// Define a function to insert a new point into a kd-tree
KDNode* 
kd_insert(KDNode* root, Point point, int depth) {
    if (root == NULL) {
        return new_kd_node(point);
    }

    // Choose the axis based on depth
    int axis = depth % 2;

    // Do not insert duplicate tuples
    if (root->point.x == point.x && root->point.y == point.y)
        return root;

    // Compare the point to the current node and recursively insert it
    if (axis == 0) {
        if (point.x < root->point.x) {
            root->left = kd_insert(root->left, point, depth+1);
        } else {
            root->right = kd_insert(root->right, point, depth+1);
        }
    } else {
        if (point.y < root->point.y) {
            root->left = kd_insert(root->left, point, depth+1);
        } else {
            root->right = kd_insert(root->right, point, depth+1);
        }
    }

    return root;
}

// Define a function to search for a point in a kd-tree
KDNode* 
kd_search(KDNode* root, Point point, int depth) {
    if (root == NULL || (root->point.x == point.x && root->point.y == point.y)) {
        return root;
    }

    // Choose the axis based on depth
    int axis = depth % 2;

    // Compare the point to the current node and recursively search for it
    if (axis == 0) {
        if (point.x < root->point.x) {
            return kd_search(root->left, point, depth+1);
        } else {
            return kd_search(root->right, point, depth+1);
        }
    } else {
        if (point.y < root->point.y) {
            return kd_search(root->left, point, depth+1);
        } else {
            return kd_search(root->right, point, depth+1);
        }
    }
}

// Define a function to search for a point in a kd-tree
inline bool
kd_exists(KDNode* root, Point point, int depth) {
    return kd_search(root, point, depth) != NULL;
}

void
kd_depth(KDNode* root, int depth, 
    int* max_depth, int* sum_depth, int *sum_nodes)
{
    int mdl = depth, mdr = depth;
    int sdl = 0, sdr = 0;
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