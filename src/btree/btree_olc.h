#ifndef __BtreeOLC_H__
#define __BtreeOLC_H__


enum class PageType : uint8_t { BTreeInner=1, BTreeLeaf=2 };

// static const uint64_t PAGESIZE=4*1024;

struct NodeBase{
  PageType type;
  uint16_t count;
};

struct BTreeLeafBase : public NodeBase {
  static const PageType typeMarker=PageType::BTreeLeaf;
};

template<class T, class Payload>
struct BTreeLeaf : public BTreeLeafBase {
  struct Entry {
    T key;
    Payload payload;
  };

   static const uint64_t maxEntries=(PAGESIZE-sizeof(NodeBase))/(sizeof(T)+sizeof(Payload));

   T keys[maxEntries];
   Payload payloads[maxEntries];
   BTreeLeaf* next;

   BTreeLeaf() {
      count=0;
      type=typeMarker;
      next = nullptr;
   }

   bool isFull() { return count==maxEntries; };

   unsigned lowerBound(T k) {
      unsigned lower=0;
      unsigned upper=count;
      do {
         unsigned mid=((upper-lower)/2)+lower;
         if (k<keys[mid]) {
            upper=mid;
         } else if (k>keys[mid]) {
            lower=mid+1;
         } else {
            return mid;
         }
      } while (lower<upper);
      return lower;
   }

   unsigned lowerBoundBF(T k) {
      auto base=keys;
      unsigned n=count;
      while (n>1) {
         const unsigned half=n/2;
         base=(base[half]<k)?(base+half):base;
         n-=half;
      }
      return (*base<k)+base-keys;
   }

  void insert(T k,Payload p) {
    assert(count<maxEntries);
    if (count) {
      unsigned pos=lowerBound(k);
      if ((pos<count) && (keys[pos]==k)) {
	// Upsert
	payloads[pos] = p;
	return;
      }
      memmove(keys+pos+1,keys+pos,sizeof(T)*(count-pos));
      memmove(payloads+pos+1,payloads+pos,sizeof(Payload)*(count-pos));
      keys[pos]=k;
      payloads[pos]=p;
    } else {
      keys[0]=k;
      payloads[0]=p;
    }
    count++;
  }

   BTreeLeaf* split(T& sep) {
      BTreeLeaf* newLeaf = new BTreeLeaf();
      newLeaf->count = count-(count/2);
      count = count-newLeaf->count;
      memcpy(newLeaf->keys, keys+count, sizeof(T)*newLeaf->count);
      memcpy(newLeaf->payloads, payloads+count, sizeof(Payload)*newLeaf->count);
      sep = keys[count-1];
      newLeaf->next = next;
      next = newLeaf;
      return newLeaf;
   }
};

struct BTreeInnerBase : public NodeBase {
  static const PageType typeMarker=PageType::BTreeInner;
};

template<class Key>
struct BTreeInner : public BTreeInnerBase {
   static const uint64_t maxEntries=(PAGESIZE-sizeof(NodeBase))/(sizeof(Key)+sizeof(NodeBase*));
   NodeBase* children[maxEntries];
   Key keys[maxEntries];

   BTreeInner() {
      count=0;
      type=typeMarker;
   }

   bool isFull() { return count==(maxEntries-1); };

   unsigned lowerBoundBF(Key k) {
      auto base=keys;
      unsigned n=count;
      while (n>1) {
         const unsigned half=n/2;
         base=(base[half]<k)?(base+half):base;
         n-=half;
      }
      return (*base<k)+base-keys;
   }

   unsigned lowerBound(Key k) {
      unsigned lower=0;
      unsigned upper=count;
      do {
         unsigned mid=((upper-lower)/2)+lower;
         if (k<keys[mid]) {
            upper=mid;
         } else if (k>keys[mid]) {
            lower=mid+1;
         } else {
            return mid;
         }
      } while (lower<upper);
      return lower;
   }

   BTreeInner* split(Key& sep) {
      BTreeInner* newInner=new BTreeInner();
      newInner->count=count-(count/2);
      count=count-newInner->count-1;
      sep=keys[count];
      memcpy(newInner->keys,keys+count+1,sizeof(Key)*(newInner->count+1));
      memcpy(newInner->children,children+count+1,sizeof(NodeBase*)*(newInner->count+1));
      return newInner;
   }

   void insert(Key k,NodeBase* child) {
      assert(count<maxEntries-1);
      unsigned pos=lowerBound(k);
      memmove(keys+pos+1,keys+pos,sizeof(Key)*(count-pos+1));
      memmove(children+pos+1,children+pos,sizeof(NodeBase*)*(count-pos+1));
      keys[pos]=k;
      children[pos]=child;
      std::swap(children[pos],children[pos+1]);
      count++;
   }

};


template<class Key,class Value>
struct BTree {
  NodeBase* root;

  BTree() {
    root = new BTreeLeaf<Key,Value>();
  }

  ~BTree() {
    delete root;
  }

  void makeRoot(Key k,NodeBase* leftChild,NodeBase* rightChild) {
    auto inner = new BTreeInner<Key>();
    inner->count = 1;
    inner->keys[0] = k;
    inner->children[0] = leftChild;
    inner->children[1] = rightChild;
    root = inner;
   }

  void insert(Key k, Value v) {
restart:
    // Current node
    NodeBase* node = root;

    // Parent of current node
    BTreeInner<Key>* parent = nullptr;

    while (node->type==PageType::BTreeInner) {
      auto inner = static_cast<BTreeInner<Key>*>(node);

      // Split eagerly if full
      if (inner->isFull()) {
        // Split
        Key sep;
        BTreeInner<Key>* newInner = inner->split(sep);
        if (parent)
          parent->insert(sep,newInner);
        else
          makeRoot(sep,inner,newInner);
        goto restart;
      }

      parent = inner;
      node = inner->children[inner->lowerBound(k)];
    }

    auto leaf = static_cast<BTreeLeaf<Key,Value>*>(node);

    // Split leaf if full
    if (leaf->count==leaf->maxEntries) {
      // Split
      Key sep;
      BTreeLeaf<Key,Value>* newLeaf = leaf->split(sep);
      if (parent)
	      parent->insert(sep, newLeaf);
      else
	      makeRoot(sep, leaf, newLeaf);
      goto restart;
    } else {
      leaf->insert(k, v);
      return; // success
    }
  }

  Value at(Key k) {
    NodeBase* node = root;
    while (node->type==PageType::BTreeInner) {
      auto inner = static_cast<BTreeInner<Key>*>(node);
      node = inner->children[inner->lowerBound(k)];
    }
    BTreeLeaf<Key,Value>* leaf = static_cast<BTreeLeaf<Key,Value>*>(node);
    unsigned pos = leaf->lowerBound(k);
    return leaf->payloads[pos];
  }

  bool exists(Key k) {
    NodeBase* node = root;
    while (node->type==PageType::BTreeInner) {
      auto inner = static_cast<BTreeInner<Key>*>(node);

      node = inner->children[inner->lowerBound(k)];
    }

    BTreeLeaf<Key,Value>* leaf = static_cast<BTreeLeaf<Key,Value>*>(node);
    unsigned pos = leaf->lowerBound(k);
    if ((pos<leaf->count) && (leaf->keys[pos]==k)) {
      return true;
    }
    return false;
  }

  uint64_t scan(Key k, int range, Value* output) {
    NodeBase* node = root;
    while (node->type==PageType::BTreeInner) {
      auto inner = static_cast<BTreeInner<Key>*>(node);

      node = inner->children[inner->lowerBound(k)];
    }

    BTreeLeaf<Key,Value>* leaf = static_cast<BTreeLeaf<Key,Value>*>(node);
    unsigned pos = leaf->lowerBound(k);
    int count = 0;
    for (unsigned i=pos; i<leaf->count; i++) {
      if (count==range)
        break;
      output[count++] = leaf->payloads[i];
    }

    return count;
  }

  std::vector<std::pair<Key, Value>> rangeQuery(Key min_key, Key max_key) {
    std::vector<std::pair<Key, Value>> result;
    NodeBase* node = root;
    while (node->type==PageType::BTreeInner) {
      auto inner = static_cast<BTreeInner<Key>*>(node);
      node = inner->children[inner->lowerBound(min_key)];
    }

    BTreeLeaf<Key,Value>* leaf = static_cast<BTreeLeaf<Key,Value>*>(node);
    bool reached_max_key = false;
    while (leaf) {
      unsigned pos = leaf->lowerBound(min_key);
      for (int i = pos; i < leaf->count; ++i) {
        if (leaf->keys[i] > max_key) {
          reached_max_key = true;
          break;
        }
        result.push_back(std::make_pair(leaf->keys[i], leaf->payloads[i]));
      }
      if (reached_max_key) {
        break;
      }
      leaf = leaf->next;
    }
    return result;
  }

  void bulk_load(std::pair<Key, uint64_t>* entries, int num_entries) {
    for (int i = 0; i < num_entries; ++i) {
      insert(entries[i].first, entries[i].second);
    }
  }

  size_t index_size() const {
    std::stack<NodeBase*> s;
    s.push(root);
    size_t size = 0;
    while (!s.empty()) {
      NodeBase* node = s.top();
      s.pop();
      size += sizeof(*node);
      if (node->type == PageType::BTreeInner) {
        auto inner = static_cast<BTreeInner<Key>*>(node);
        for (int i = 0; i < inner->count; ++i) {
          s.push(inner->children[i]);
        }
      }
    }
    return size;
  }

};

#endif // __BtreeOLC_H__