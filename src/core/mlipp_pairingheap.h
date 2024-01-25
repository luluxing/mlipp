#ifndef __PAIRINGHEAP_H_
#define __PAIRINGHEAP_H_

#include "point.h"

typedef struct MPNode {
    struct MPNode *first_child;
    struct MPNode *next_sibling;
    struct MPNode *prev_or_parent;
    double distance;
    double dx;
    double dy;
    int    entry_type;
    void  *node;
    int    pos;
} MPNode;

typedef struct MPHeap {
    MPNode *root;
} MPHeap;


#define MP_LEFT  1
#define MP_RIGHT 2
#define MP_ANY   3

#define pairingheap_is_empty(heap) ((heap->root) == NULL)

/*
 * pairingheap_allocate
 *
 * Returns a pointer to a newly-allocated heap
 */
MPNode *
pairingheap_allocate(int entry_type, void *node, int pos,
    double dx, double dy, double distance) {
    MPNode *root = (MPNode*)malloc(sizeof(MPNode));
    root->first_child = NULL;
    root->next_sibling = NULL;
    root->prev_or_parent = NULL;
    root->distance = distance;
    root->dx = dx;
    root->dy = dy;
    root->entry_type = entry_type;
    root->node = node;
    root->pos = pos;
    return root;
}

/*
 * A helper function to merge two subheaps into one.
 *
 * The subheap with smaller value is put as a child of the other one (assuming
 * a max-heap).
 *
 * The next_sibling and prev_or_parent pointers of the input nodes are
 * ignored. On return, the returned node's next_sibling and prev_or_parent
 * pointers are garbage.
 */
static MPNode *
merge(MPNode *a, MPNode *b)
{
    if (a == NULL)
        return b;
    if (b == NULL)
        return a;

    /* swap 'a' and 'b' so that 'a' is the one with larger value */
    if (a->distance > b->distance)
    {
        MPNode *tmp;

        tmp = a;
        a = b;
        b = tmp;
    }

    /* and put 'b' as a child of 'a' */
    if (a->first_child)
        a->first_child->prev_or_parent = b;
    b->prev_or_parent = a;
    b->next_sibling = a->first_child;
    a->first_child = b;
    return a;
}

/*
 * pairingheap_insert
 *
 * Adds the given node to the heap in O(1) time.
 */
void 
pairingheap_add(MPHeap *heap, MPNode *node) {
    node->first_child = NULL;
    /* Link the new node as a new tree */
    heap->root = merge(heap->root, node);
    heap->root->prev_or_parent = NULL;
    heap->root->next_sibling = NULL;
}

/*
 * Merge a list of subheaps into a single heap.
 *
 * This implements the basic two-pass merging strategy, first forming pairs
 * from left to right, and then merging the pairs.
 */
static MPNode *
merge_children(MPNode *children)
{
    MPNode *curr, *next;
    MPNode *pairs;
    MPNode *newroot;

    if (children == NULL || children->next_sibling == NULL)
        return children;

    /* Walk the subheaps from left to right, merging in pairs */
    next = children;
    pairs = NULL;
    for (;;)
    {
        curr = next;

        if (curr == NULL)
            break;

        if (curr->next_sibling == NULL)
        {
            /* last odd node at the end of list */
            curr->next_sibling = pairs;
            pairs = curr;
            break;
        }

        next = curr->next_sibling->next_sibling;

        /* merge this and the next subheap, and add to 'pairs' list. */

        curr = merge(curr, curr->next_sibling);
        curr->next_sibling = pairs;
        pairs = curr;
    }

    /*
     * Merge all the pairs together to form a single heap.
     */
    newroot = pairs;
    next = pairs->next_sibling;
    while (next)
    {
        curr = next;
        next = curr->next_sibling;

        newroot = merge(newroot, curr);
    }

    return newroot;
}

/*
 * pairingheap_remove_first
 *
 * Removes the first (root, topmost) node in the heap and returns a pointer to
 * it after rebalancing the heap. The caller must ensure that this routine is
 * not used on an empty heap. O(log n) amortized.
 */
bool
pairingheap_remove_first(MPHeap *heap, MPNode **node)
{
    if (pairingheap_is_empty(heap))
        return false;

    /* Remove the root, and form a new heap of its children. */
    MPNode *children = heap->root->first_child;
    *node = heap->root;

    heap->root = merge_children(children);
    if (!pairingheap_is_empty(heap))
    {
        heap->root->prev_or_parent = NULL;
        heap->root->next_sibling = NULL;
    }
    return true;
}

#endif // __PAIRINGHEAP_H_