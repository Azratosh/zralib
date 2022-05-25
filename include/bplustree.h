#ifndef BPLUSTREE_H
#define BPLUSTREE_H

/*
This entire header stems from a university project and needs to be refactored.

Read at your own peril.
*/

#include <algorithm>
#include <cstddef>
#include <exception>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>

#ifdef DEBUG
#define DEBUG_STDERR(x) (std::cerr << x << '\n')
#define DEBUG_STDOUT(x) (std::cout << x << '\n')
#else
#define DEBUG_STDERR(x)                                                        \
  do {                                                                         \
  } while (0)
#define DEBUG_STDOUT(x)                                                        \
  do {                                                                         \
  } while (0)
#endif

namespace zra {

/*
Default properties for each tree that is created.
Note that root nodes are an exception to these.
*/
template <typename Key, size_t N> struct bplustree_properties {
  const static size_t min_store{N}; // ! root: at least 1 key
  const static size_t max_store{N * 2};
  const static size_t min_child_nodes{N + 1}; // ! root: at least 2 child nodes
  const static size_t max_child_nodes{N * 2 + 1};
};

/*
Implementation of a B+-tree for Algorithms and Data Structures (ADS)
at the University of Vienna, made in the winter semester of 2021.
*/
template <typename Key, size_t N = 10, size_t MaxCacheSize = 16,
          typename Properties = bplustree_properties<Key, N>>
class bplustree {
  // specification, type aliases
public:
  class ConstIterator;
  // class ConstReverseIterator;

  using value_type = Key;
  using reference = value_type &;
  using const_reference = const value_type &;
  using pointer = value_type *;
  using const_pointer = const value_type *;

  using key_type = Key;

  using size_type = size_t;
  using difference_type = std::ptrdiff_t;
  using key_compare = std::less<key_type>;
  using key_equal_to = std::equal_to<key_type>;

  using const_iterator = bplustree::ConstIterator;

  // ! This is done so as to match the specification on CEWebs
  using iterator = const_iterator;

  // ? These should be implemented at some point as well?
  // using const_reverse_iterator = ADS_set::ConstReverseIterator;
  // using reverse_iterator = ADS_set::ReverseIterator;

  /*
  The nodes of the B+-Tree.

  IndexNode and LeafNode both inherit from Node, which is an abstract struct.
  The Node struct itself should never be instantiated by itself and is only
  to be used for parametric polymorphism using Node::is_leaf();
  LeafNodes namely always have a `layer` of `0`.

  TODO: insert, erase for Nodes
  TODO: fix bplustree::delete() being slow as fuck in certain cases
  */
private:
  struct Node {
    // Nodes with layer 0 must always be leaf nodes
    size_type layer;
    size_type size;

    Node(const size_type &layer) : layer{layer}, size{0} {}
    virtual ~Node() = default;

    const key_type at(const size_type &) const { return key_type(); }

    virtual void insert(const size_type &, const key_type &) {}
    virtual void insert(const size_type &, const key_type &, Node *) {}

    virtual void reset() { this->size = 0; }

    virtual bool is_below_min() const {
      return this->size < Properties::min_store;
    }

    virtual bool is_below_or_at_min() const {
      return this->size <= Properties::min_store;
    }

    virtual bool is_above_min() const {
      return Properties::min_store < this->size &&
             this->size <= Properties::max_store;
    }

    virtual bool is_at_max() const {
      return this->size == Properties::max_store;
    }

    virtual bool is_leaf() const { return layer == 0; }

    virtual std::ostream &print(std::ostream &o) const {
      return o << "Node {l: " << this->layer << ", s: " << this->size << "}[]";
    }

    // Operators

    inline friend std::ostream &operator<<(std::ostream &o, const Node *&node) {
// Allow printing nullptrs in debug mode
#ifdef DEBUG
      if (node == nullptr) {
        return o << "NULLPTR";
      }
#endif

      if (node->is_leaf()) {
        return dynamic_cast<const LeafNode *>(node)->print(o);
      } else {
        return dynamic_cast<const IndexNode *>(node)->print(o);
      }
    }
  };

  struct IndexNode : public Node {
    // Note: Keys equal values, keep that in mind

    key_type keys[Properties::max_store];
    Node *children[Properties::max_child_nodes];

    IndexNode(const size_type &layer = 1) : Node(layer) {}

    const key_type &at(const size_type &pos) const {
      if (pos < this->size) {
        return this->keys[pos];
      }

      throw std::runtime_error("Key does not exist");
    }

    void insert(const size_type &pos, const key_type &key,
                Node *const &child_node) {
      DEBUG_STDOUT("Inserting " << key << " and " << child_node
                                << " at position " << pos);

      if (child_node == nullptr) {
        throw std::runtime_error("Received nullptr for key at position " +
                                 std::to_string(pos));
      }

      if (pos > this->size) {
        throw std::runtime_error(
            "Position out of bounds: size = " + std::to_string(this->size) +
            ", pos: " + std::to_string(pos));
      }
      // ! Special case: extra child node shifted separately first
      this->children[this->size + 1] = this->children[this->size];

      // ! NOTE: Values are *always* shifted because IndexNode
      // !       gets Child Node assigned at pos 0 on creation
      // Shift all keys and children from insert pos to the right
      for (size_type i{this->size}; i > pos; --i) {
        this->keys[i] = this->keys[i - 1];
        this->children[i] = this->children[i - 1];
      }

      this->keys[pos] = key;
      this->children[pos + 1] = child_node;
      this->size++;

      DEBUG_STDOUT("Inserted key " << key << ", child " << child_node
                                   << " at position " << pos);
    }

    void reset() { Node::reset(); }

    std::ostream &print(std::ostream &o) const {
      o << "IndexNode {l: " << this->layer << ", s: " << this->size
        << ", c: " << (this->size > 0 ? this->size + 1 : 0) << "}[";

      if (this->size == 1) {
        o << this->keys[0] << "](0, 1)";

      } else if (this->size > 1) {
        o << this->keys[0];
        for (size_type i{1}; i < this->size; ++i) {
          o << ", " << this->keys[i];
        }

        o << "](0";
        for (size_type i{1}; i < this->size + 1; ++i) {
          o << ", " << i;
        }
        o << ")";

      } else {
        o << "]()";
      }

      return o;
    }
  };

  struct LeafNode : public Node {
    // Same as linked list
    LeafNode *left_leaf;
    LeafNode *right_leaf;
    key_type values[Properties::max_store];

    LeafNode() : Node(0), left_leaf{nullptr}, right_leaf{nullptr} {}

    const key_type &at(const size_type &pos) const {
      if (pos < this->size) {
        return this->values[pos];
      }

      throw std::runtime_error("Value does not exist");
    }

    void insert(const size_type &pos, const key_type &key) {
      if (pos > this->size) {
        throw std::runtime_error(
            "Position out of bounds: size = " + std::to_string(this->size) +
            ", pos: " + std::to_string(pos));
      }

      // size_type last_element{this->size};
      // while (last_element != insert_position) {
      //   last_element--;
      //   this->values[last_element] = this->values[last_element - 1];
      // }

      // Shift all values from insert pos to the right
      for (size_type i{this->size}; i > pos; --i) {
        this->values[i] = this->values[i - 1];
      }

      this->values[pos] = key;
      this->size++;

      DEBUG_STDOUT("Inserted value " << key << " at position " << pos);
    }

    void reset() {
      Node::reset();
      this->right_leaf = nullptr;
      this->left_leaf = nullptr;
    }

    std::ostream &print(std::ostream &o) const {
      o << "LeafNode {l: " << this->layer << ", s: " << this->size << "}[";

      if (this->size > 0) {
        o << this->values[0];
        for (size_t i{1}; i < this->size; ++i) {
          o << ", " << this->values[i];
        }
      }

      return o << "]";
    }
  };

  /*
  Keeps nodes locked up until they're ready to be used.

  In other words, nodes that would be deallocated are put
  back into the dungeon; as long as nodes are available,
  they're freed and reused instead of allocated.

  * Basically a stack that recycles nodes.
  */
  template <typename node_type> class NodeDungeon {
  private:
    size_type _size{};
    node_type *_node_stack[MaxCacheSize];

    void _allocate() {
      node_type *_new_node_ptr = new node_type();
      this->_node_stack[this->_size] = _new_node_ptr;
      this->_size++;
    }

  public:
    size_type size() const { return this->_size; }

    ~NodeDungeon() { this->free(); }

    /*
    Returns the last node pointer in the stack. If no node
    is available, a new one is created instead.
    */
    node_type *pop() {
      if (this->_size == 0) {
        return new node_type();
      }

      return this->_node_stack[--this->_size];
    }

    /*
    Adds a node pointer to the stack. If the stack
    is full, the node is deleted instead.
    */
    void push(node_type *const &node) {
      if (_size < MaxCacheSize) {
        node->reset();

        this->_node_stack[_size] = node;
        _size++;

      } else {
        delete node;
      }
    }

    /*
    Deletes all nodes on the stack, setting the size
    of the stack to `0`.
    */
    void free() {
      for (size_type i{this->_size}; i-- > 0;) {
        delete this->_node_stack[i];
      }
      this->_size = 0;
    }

    /*
    Pre-allocate a bunch of nodes in memory.
    */
    void allocate(size_type n) {
      size_type max_allocatable{(MaxCacheSize - this->_size) - 1};
      if (n > max_allocatable) {
        n = max_allocatable;
      }

      if (n == 0) {
        return;
      }

      for (size_type i{n}; i-- > 0;) {
        this->_allocate();
      }
    }
  };

  /*
  The B+-Tree's instance variables.

  It's only necessary to keep track of the root node and outermost leaf nodes.
  The root node is used for traversion, whereas the two leaf nodes are used
  for iteration.
  */
private:
  Node *_root_node;
  LeafNode *_leftmost_leaf;
  LeafNode *_rightmost_leaf;

  NodeDungeon<LeafNode> *_l_node_cache;
  NodeDungeon<IndexNode> *_i_node_cache;

  size_type _size;
  size_type _leaf_nodes;
  size_type _index_nodes;

public:
  // * * * Constructors * * *
  bplustree() : _size{0}, _leaf_nodes{1}, _index_nodes{0} {
    this->_l_node_cache = new NodeDungeon<LeafNode>();
    this->_i_node_cache = new NodeDungeon<IndexNode>();

    this->_root_node = this->_leftmost_leaf = this->_rightmost_leaf =
        this->_l_node_cache->pop();
  }

  bplustree(std::initializer_list<key_type> ilist) : bplustree() {
    this->insert(ilist);
  }

  template <typename InputIt>
  bplustree(InputIt first, InputIt last) : bplustree() {
    this->insert(first, last);
  }

  bplustree(const bplustree &other) : bplustree() {
    for (const auto &element : other) {
      this->insert(element);
    }
  }

  ~bplustree() {
    DEBUG_STDOUT("Destructing tree");

    this->clear();

    DEBUG_STDOUT("Deleting root node");

    delete this->_root_node;

    delete this->_l_node_cache;
    delete this->_i_node_cache;

    DEBUG_STDOUT("Deleted root node; destructed tree");
  }

  // * * * Operators * * *

  bplustree &operator=(const bplustree &other) {
    if (this != &other) {
      this->clear();
      for (const auto &element : other) {
        this->insert(element);
      }
    }
    return *this;
  }

  bplustree &operator=(std::initializer_list<key_type> ilist) {
    this->clear();
    this->insert(ilist);
    return *this;
  }

  // * * * Methods * * *

  size_type size() const { return this->_size; }

  size_type index_nodes() const { return this->_index_nodes; }

  size_type leaf_nodes() const { return this->_leaf_nodes; }

  size_type nodes() const { return this->_index_nodes + this->_leaf_nodes; }

  bool empty() const { return this->_size == 0; }

  void insert(std::initializer_list<key_type> ilist) {
    this->allocate_in_cache(ilist.size());

    for (const auto &element : ilist) {
      this->_init_insert(element);
    }
  }

  std::pair<iterator, bool> insert(const key_type &key) {
    return this->_init_insert(key);
  }

  template <typename InputIt> void insert(InputIt first, InputIt last) {
    std::for_each(first, last,
                  [&](const auto &element) { this->_init_insert(element); });
  }

  /*
  Allocates enough extra nodes in the cache corresponding to
  the number of keys `key_count` to hold.

  If there are still enough nodes available, no new nodes are allocated.
  */
  void allocate_in_cache(size_type key_count) {
    if (key_count < Properties::max_store) {
      return;
    }

    size_type pre_alloc_leaves{
        key_count / (Properties::min_store + (Properties::min_store / 2))};

    if (pre_alloc_leaves > 0 &&
        this->_l_node_cache->size() < pre_alloc_leaves) {
      this->_l_node_cache->allocate(pre_alloc_leaves -
                                    this->_l_node_cache->size());
    }

    size_type pre_alloc_indices{
        pre_alloc_leaves /
        (Properties::min_child_nodes + (Properties::min_child_nodes / 2))};

    if (pre_alloc_indices > 0 &&
        this->_i_node_cache->size() < pre_alloc_indices) {
      this->_i_node_cache->allocate(pre_alloc_indices -
                                    this->_i_node_cache->size());
    }
  }

  /*
  Yeets all the nodes recursively and resets
  root, leftmost leaf, rightmost leaf, and all stats.
  */
  void clear() {
    if (this->_root_node != nullptr) {
      DEBUG_STDOUT("Clearing ...");

      this->_clear(this->_root_node);

      DEBUG_STDOUT("Deleting root node");
      if (this->_root_node->is_leaf()) {
        this->_l_node_cache->push(dynamic_cast<LeafNode *>(this->_root_node));
      } else {
        this->_i_node_cache->push(dynamic_cast<IndexNode *>(this->_root_node));
      }
      DEBUG_STDOUT("Deleted root node");

      DEBUG_STDOUT("Allocating new root node");
      this->_root_node = this->_leftmost_leaf = this->_rightmost_leaf =
          this->_l_node_cache->pop();
      DEBUG_STDOUT("Allocated new root node");

      this->_size = 0;
      this->_index_nodes = 0;
      this->_leaf_nodes = 1;
    }
  }

  /*
  Returns number of deleted keys (0 or 1).
  */
  size_type erase(const key_type &key) {
    auto r{this->_erase_init(key)};

#ifdef DEBUG
    this->dump_detailed();
#endif

    return r;
  }

  void swap(bplustree &other) {
    std::swap(this->_size, other._size);
    std::swap(this->_index_nodes, other._index_nodes);
    std::swap(this->_leaf_nodes, other._leaf_nodes);

    std::swap(this->_root_node, other._root_node);
    std::swap(this->_leftmost_leaf, other._leftmost_leaf);
    std::swap(this->_rightmost_leaf, other._rightmost_leaf);

    std::swap(this->_i_node_cache, other._i_node_cache);
    std::swap(this->_l_node_cache, other._l_node_cache);
  }

  iterator find(const key_type &key) const {
#ifdef DEBUG
    this->dump_detailed();
#endif

    DEBUG_STDOUT("FIND: Looking for key " << key);

    if (this->_size == 0) {
      DEBUG_STDOUT("FIND ERROR: Early abort: size == 0");
      return this->end();
    }

    Node *current_node = this->_root_node;
    while (!current_node->is_leaf()) {
      IndexNode *i_node{dynamic_cast<IndexNode *>(current_node)};

      DEBUG_STDOUT("-- Depth: " << this->_root_node->layer - current_node->layer
                                << " - IndexNode: " << i_node);

      size_type i_position{this->get_pos(i_node, key)};
      if (i_position < i_node->size &&
          key_equal_to{}(key, i_node->at(i_position))) {
        ++i_position;
      }

      DEBUG_STDOUT("-- Encountered Node at [" << i_position << "]");

      current_node = i_node->children[i_position];
    }

    LeafNode *l_node{dynamic_cast<LeafNode *>(current_node)};
    size_type l_position{this->get_pos(l_node, key)};

    DEBUG_STDOUT("-- Depth: " << this->_root_node->layer - current_node->layer
                              << " - LeafNode: " << l_node);

    if (l_position < l_node->size &&
        key_equal_to{}(key, l_node->values[l_position])) {
      DEBUG_STDOUT("FIND SUCCESS: " << key << " at position " << l_position
                                    << ": " << l_node);

      return iterator(l_node, l_position);
    } else {
      DEBUG_STDOUT("FIND ERROR: Couldn't find key " << key);

      return this->end();
    }
  }

  /*
  This is a *really* lazy implementation, but hey,
  I'm not gonna write the same thing twice
  */
  size_type count(const key_type &key) const {
    if (this->_size == 0) {
      return 0;
    }

    return static_cast<size_type>(this->find(key) != this->end());
  }

  iterator begin() const { return iterator(this->_leftmost_leaf); }

  iterator end() const {
    return iterator(this->_rightmost_leaf, this->_rightmost_leaf->size);
  }

  // debugging
  void dump(std::ostream &o = std::cerr) const {
    o << std::endl;
    o << "\n= = = = = = = = = = = = = = = = = = = =\n";
    o << "DUMP - Dumping tree:\n";
    o << "Size:      " << this->_size << "\n";

    o << "Elements: {";
    for (const auto &element : *this) {
      o << " " << element;
    }
    o << " }\n";

    o << "Leaves:\n";
    LeafNode *_leaf = this->_leftmost_leaf;
    do {
      o << "\t" << _leaf << "\n";
      _leaf = _leaf->right_leaf;
    } while (_leaf != nullptr);

#ifdef DEBUG
    this->dump_detailed(o);
#endif

    o << "DUMP - Ending dump\n";
    o << "\n= = = = = = = = = = = = = = = = = = = =\n";
  }

  // * * * Private Methods * * *

private:
  void _clear(Node *const &current_node) {
    if (!current_node->is_leaf()) {
      IndexNode *i_node{dynamic_cast<IndexNode *>(current_node)};
      for (size_type i{0}; i < i_node->size + 1; ++i) {
        DEBUG_STDOUT("-- Depth: " << (this->_root_node->layer - i_node->layer)
                                  << " | Descending into: "
                                  << i_node->children[i]);

        this->_clear(i_node->children[i]);

        DEBUG_STDOUT("-- Depth: " << (this->_root_node->layer - i_node->layer)
                                  << " | Deleting Node: "
                                  << i_node->children[i]);

        if (i_node->children[i]->is_leaf()) {
          this->_l_node_cache->push(
              dynamic_cast<LeafNode *>(i_node->children[i]));
        } else {
          this->_i_node_cache->push(
              dynamic_cast<IndexNode *>(i_node->children[i]));
        }
      }
    }
  }

  void dump_detailed(std::ostream &o = std::cerr,
                     Node *const &passed_node = nullptr,
                     const size_type &current_iteration = 0) const {
    const Node *current_node{passed_node == nullptr ? this->_root_node
                                                    : passed_node};
    size_type depth{this->_root_node->layer - current_node->layer};

    if (current_iteration == 0) {
      o << "\n=== DETAILED DUMP ===\n";
    }

    for (size_type i{0}; i < depth; ++i) {
      o << "--";
    }
    o << " [" << current_iteration << "] " << current_node << "\n";

    if (current_node->size != 0 && !current_node->is_leaf()) {
      const IndexNode *index_node{
          dynamic_cast<const IndexNode *>(current_node)};
      if (index_node->children[0]->is_leaf()) {

        for (size_type i{0}; i < index_node->size + 1; ++i) {
          o << "  ";
          for (size_type j{0}; j < depth; ++j) {
            o << "--";
          }
          o << "[" << i << "] " << index_node->children[i] << "\n";
        }
      } else {
        for (size_type i{0}; i < index_node->size + 1; ++i) {
          this->dump_detailed(o, index_node->children[i],
                              current_iteration + 1);
        }
      }
    }

    if (current_iteration == 0) {
      o << "\n=== DETAILED DUMP END ===\n";
    }
  }

  /*
  Gets the position of a key in the given node,
  which may be used to either traverse nodes, or find
  a valid position to insert or erase the given key.

  The maximum possible value returned is equal to the size of the given node.
  This can be leveraged to traverse `IndexNode`s.

  ! Does not check for equality.

  For insert, find, erase, etc. equality must still be checked.

  * This could be optimized with binary search, but only for large N
  */
  template <typename node_pointer>
  size_type get_pos(const node_pointer *node, const key_type &key) const {
    DEBUG_STDOUT("GET_POS: Getting position for key " << key << " in " << node);

    size_type pos{0};

    while (pos < node->size && key_compare{}(node->at(pos), key)) {
      ++pos;
    }

    DEBUG_STDOUT("GET_POS: Position for key: " << pos);

    return pos;
  }

  // * * * Insertion * * *

  std::pair<iterator, bool> _init_insert(const key_type &key) {
    Node *propagated_node = nullptr;
    key_type propagated_key = key_type();

    DEBUG_STDOUT("Recursive insertion for: " << key);

    std::pair<iterator, bool> result{this->_insert(
        this->_root_node, key, &propagated_node, &propagated_key)};

    // Handling root overflow
    // This is the only kind of overflow that isn't handled in _insert
    // as it is a special case of overflowing
    if (propagated_node != nullptr) {
      // ! Layer is incremented since this is a root split

      DEBUG_STDOUT("Creating new root node");
      IndexNode *new_root_node = this->_i_node_cache->pop();
      new_root_node->layer = this->_root_node->layer + 1;

      new_root_node->children[0] = this->_root_node;
      new_root_node->insert(0, propagated_key, propagated_node);

      this->_root_node = new_root_node;

      DEBUG_STDOUT("Created new root node: " << new_root_node);
      DEBUG_STDOUT("- Child 0: " << new_root_node->children[0]);
      DEBUG_STDOUT("- Child 1: " << new_root_node->children[1]);
    }

    DEBUG_STDOUT("Insertion result: " << *result.first << ": "
                                      << result.second);

#ifdef DEBUG
    this->dump_detailed();
#endif
    return result;
  }

  /*
  Recursively moves down the tree and inserts the key where possible.
  If the key was inserted successfully, the `bool` provided with `iterator`
  will be `true`. If the key already exists, it will be `false` instead.
  */
  std::pair<iterator, bool> _insert(Node *const node, const key_type &key,
                                    Node **const propagated_node_ptr,
                                    key_type *propagated_key_ptr) {
    /*
    The entire tree is traversed down to the key's corresponding leaf.
    If the leaf node is split, the split is trickled upwards again,
    which may trigger multiple index node splits.
    */

    DEBUG_STDOUT("-- Depth: " << (this->_root_node->layer - node->layer)
                              << " | Inserting " << key);

    if (node->is_leaf()) {
      // 0: Find position to insert key
      LeafNode *l_node{dynamic_cast<LeafNode *>(node)};
      size_type insert_position{this->get_pos(l_node, key)};

      DEBUG_STDOUT("Current Node: " << l_node
                                    << " Key Position: " << insert_position);

      // 1: Return iterator + false if key already exists
      if (insert_position < l_node->size &&
          key_equal_to{}(key, l_node->values[insert_position])) {
        return std::pair<iterator, bool>(iterator(l_node, insert_position),
                                         false);
      }

      // 2: Handle overflow
      if (l_node->is_at_max()) {
        DEBUG_STDOUT("Handling LeafNode overflow");

        // 2.1 Split the leaf and get its new right sibling
        LeafNode *new_l_node{this->split_node(l_node)};

        // 2.2 If position is in next node, adjust position
        // and reassign leaf node to save code
        if (l_node->size <= insert_position) {
          insert_position -= l_node->size;
          l_node = new_l_node;
        }

        // 2.3: Propagate the split upwards; the new node and key
        // will be inserted in the parent index node
        *propagated_node_ptr = new_l_node;
      }

      // 3: Insert the key into the node at the given position
      l_node->insert(insert_position, key);
      DEBUG_STDOUT("LeafNode after insertion: " << l_node);

      this->_size++;
      DEBUG_STDOUT("Increasing size to " << this->_size);

      // 4: If there's a new LeafNode, set its first value as propagated key
      if (*propagated_node_ptr != nullptr) {
        *propagated_key_ptr =
            dynamic_cast<LeafNode *>(*propagated_node_ptr)->values[0];
      }

      return std::pair<iterator, bool>(iterator(l_node, insert_position), true);

    } else { // node == IndexNode
      // 0: Find position to insert key
      IndexNode *i_node{dynamic_cast<IndexNode *>(node)};
      size_type insert_position{this->get_pos(i_node, key)};

      DEBUG_STDOUT("Current Node: " << i_node
                                    << " Key Position: " << insert_position);

      // ! If key at insert pos is equal to key to insert, increment pos
      if (insert_position < i_node->size &&
          key_equal_to{}(key, i_node->at(insert_position))) {
        ++insert_position;
      }

      Node *new_child_node_ptr = nullptr;
      key_type new_propagated_key = key_type();

      // 1: Traverse down recursively until leaf node is encountered
      // ! insert_position will be equal to node's size if no key was found
      // Using this behaviour makes it possible to access the last child
      // node
      std::pair<iterator, bool> result{
          _insert(i_node->children[insert_position], key, &new_child_node_ptr,
                  &new_propagated_key)};
      // - Split index node

      // 2: Catch new child node that was created on split
      if (new_child_node_ptr != nullptr) {
        // 2.1: Handle overflow of index node
        DEBUG_STDOUT("Encountered new child node");

        if (i_node->is_at_max()) {
          DEBUG_STDOUT("Handling IndexNode overflow");

          // 2.1.a: Split index node and get its new sibling as well
          // as the key for the parent index node
          auto [new_i_node, new_key] =
              this->split_node(i_node, insert_position);

          // 2.1.b: Handle insertion right at split position
          // * See split_node() for IndexNode
          if (i_node->size + 1 == insert_position &&
              i_node->size < new_i_node->size) {
            DEBUG_STDOUT("Handling insertion at split position");

            // ! New key gets inserted at the back instead of propagated key
            i_node->keys[i_node->size] = new_key;

            // First child node of new index node gets inserted at
            // end of original node
            i_node->children[i_node->size + 1] = new_i_node->children[0];

            // ... and the new child node gets inserted
            // at the front of the new node instead
            new_i_node->children[0] = new_child_node_ptr;

            // Increase size and propagate node upwards; early return
            i_node->size++;
            *propagated_node_ptr = new_i_node;
            // ! Propagated key is not changed and moves further up
            // I don't fucking know why I need to reassign it here
            *propagated_key_ptr = new_propagated_key;

            DEBUG_STDOUT(
                "Original IndexNode (split position handling): " << i_node);
            DEBUG_STDOUT("Propagated node: " << *propagated_node_ptr);
            DEBUG_STDOUT("Propagated key: " << *propagated_key_ptr);

            return result;
          }

          // 2.1.c: Reassign current index node if insert position
          // is in the newly split node to save code
          if (i_node->size < insert_position) {
            DEBUG_STDOUT("Handling insert position "
                         << insert_position
                         << " being in new node; size of i_node being "
                         << i_node->size);

            insert_position -= i_node->size + 1;
            i_node = new_i_node;

            DEBUG_STDOUT("Reassigned node, adjusted position to "
                         << insert_position);
          }

          // 2.1.d: Propagate the split upwards; the new node and key
          // will be inserted in the parent index node
          *propagated_node_ptr = new_i_node;
          *propagated_key_ptr = new_key;

          DEBUG_STDOUT("Propagated Node: " << *propagated_node_ptr
                                           << "\nPropagated Key: "
                                           << *propagated_key_ptr);
        }

        // 3: Finally, put the key and child node into
        // the index node at the given position
        i_node->insert(insert_position, new_propagated_key, new_child_node_ptr);
      }

      DEBUG_STDOUT("IndexNode after insertion: " << i_node);
      return result;
    }
  }

  /*
  Split a leaf node and return a pointer to the newly created node,
  as well as the key for the nodes' parent index node.
  */
  LeafNode *split_node(LeafNode *&l_node) {
    // Note: New LeafNodes are inserted to the
    //       "right" of the given node

    // Note: Unlike when splitting index nodes, the middle of leaf
    //       nodes doesn't need to be adjusted

    DEBUG_STDOUT("SPLIT: " << l_node);

    // 0: Find middle
    size_type middle{l_node->size / 2};

    // 1: Create new leaf
    LeafNode *new_l_node = this->_l_node_cache->pop();

    // Set size of new node to (about) half of the size of the original node
    // ! Necessary if for some reason we want to split only partially filled
    // ! nodes later
    new_l_node->size = l_node->size - middle;

    // 2: Modify leaf links
    // Right node of original becomes right node of new node
    new_l_node->right_leaf = l_node->right_leaf;

    // Handle rightmost (leftmost should never change during split)
    if (new_l_node->right_leaf == nullptr) {
      // New node becomes rightmost leaf of linked list of leaves
      this->_rightmost_leaf = new_l_node;
    } else {
      // New node becomes the left sibling of its right sibling
      new_l_node->right_leaf->left_leaf = new_l_node;
    }

    // New node becomes right sibling of the original node
    l_node->right_leaf = new_l_node;

    // Original node becomes left sibling of the new node
    new_l_node->left_leaf = l_node;

    // * At this point, the new has been correctly created
    // * and linked with its siblings

    // 3: Copy values to new leaf node from middle onwards
    for (size_type i{middle}; i < l_node->size; ++i) {
      new_l_node->values[i - middle] = l_node->values[i];
    }

    // 4: Resize original leaf node
    l_node->size = middle;

    // 6: Woo, stats
    this->_leaf_nodes++;

    DEBUG_STDOUT("Original Node: " << l_node);
    DEBUG_STDOUT("New Node:      " << new_l_node);
    DEBUG_STDOUT("SPLIT END");

    return new_l_node;
  }

  /*
  Splits an index node and returns a pointer to the newly created node,
  as well as the middle key for the nodes' parent index node.

  An additional insertion position must be provided in order to
  handle uneven splits.
  */
  std::pair<IndexNode *, key_type>
  split_node(IndexNode *&i_node, const size_type &insert_position) {
    // Note: New IndexNodes are inserted to the
    // "right" too

    DEBUG_STDOUT("SPLIT: " << i_node);

    // 0: Find middle
    size_type middle{i_node->size / 2};

    // If insert position is below or at middle,
    // the middle must be offset by -1 to account for key propagation
    // Otherwise, under- / overflows might happen
    if (insert_position <= middle) {
      DEBUG_STDOUT("Handling insertion before or at middle");
      middle--;
    }

    // 1: Create new index node
    IndexNode *new_i_node{this->_i_node_cache->pop()};
    new_i_node->layer = i_node->layer;

    // ! Offset of -1 again is necessary for key propagation
    // ! The actual middle key is put into the parent index node
    new_i_node->size = (i_node->size - middle) - 1;

    // 2: Copy every key after the middle to new index node
    for (size_type i{middle + 1}; i < i_node->size; ++i) {
      new_i_node->keys[i - middle - 1] = i_node->keys[i];
    }

    // 3: Copy every child after the middle to the new index node
    for (size_type i{middle + 1}; i < i_node->size + 1; ++i) {
      new_i_node->children[i - middle - 1] = i_node->children[i];
    }

    // 4: Get the middle key to propagate to the parent index node
    key_type middle_key{i_node->keys[middle]};

    // 5: Resize original node
    /*
    ! Important: If middle was decremented because it was uneven and the new
    ! key will be inserted into the original node, the original node's size
    ! will be new node's size - 1 at this point. After insertion, their
    ! sizes will be equal.
    */
    i_node->size = middle;

    // 6: Woo, more stats
    this->_index_nodes++;

    DEBUG_STDOUT("Original Node: " << i_node);
    DEBUG_STDOUT("New Node:      " << new_i_node);
    DEBUG_STDOUT("Middle Key:    " << middle_key);

    return std::pair<IndexNode *, key_type>(new_i_node, middle_key);
  }

  // * * * Deletion * * *

  size_type _erase_init(const key_type &key) {
    DEBUG_STDOUT("ERASING: " << key);

    if (this->_size == 0) {
      DEBUG_STDOUT("Size is 0, nothing to erase");
      return 0;
    }

    DelState deletion_state{this->_erase(key, this->_root_node)};

    DEBUG_STDOUT("Received DelState: " << deletion_state);

    if (deletion_state == DelState::SUCCESS) {
      this->_size--;
      return 1;
    } else if (deletion_state == DelState::FAILURE) {
      return 0;
    } else {
      throw std::runtime_error(
          "Passed merge to _erase_init -- how did that happen?");
    }

    DEBUG_STDOUT("ERASION END: " << key);
  }

  enum DelState : short {
    SUCCESS, // * Deletion succeeded    => nothing to do
    FAILURE, // * Deletion failed       => key not found
    MERGE,   // * Deletion caused merge => parent node must delete empty node
    NON_INIT = 100 // * Special state to make the compiler shut up
  };

  enum DelCase : short {
    // ! LeafNode: do nothing; IndexNode: Delete if only one child
    HANDLE_ROOT,

    // * Node to the right has extra data, balance to the left with current
    BALANCE_RIGHT_TO_CURRENT,

    // * Node to the left has extra data, balance to the right with current
    BALANCE_LEFT_TO_CURRENT,

    // * Both left and right node would underflow, merge left and current
    MERGE_LEFT_AND_CURRENT,

    // * Both left and right node would underflow, merge right and current
    MERGE_RIGHT_AND_CURRENT,

    // * Nothing needs to be handled
    NO_ACTION
  };

  DelCase _determine_delcase(const Node *const &current_node,
                             const Node *const &current_node_left,
                             const Node *const &current_node_right,
                             const IndexNode *const &parent_node,
                             const IndexNode *const &parent_node_left,
                             const IndexNode *const &parent_node_right) const {
    DEBUG_STDOUT("DELCASE HANDLING");

    // ! Case 0: Nothing needs to be done
    if (!current_node->is_below_min()) {
      DEBUG_STDOUT("NO_ACTION: Node is not below minimum");
      return DelCase::NO_ACTION;
    }

    // ! Case 1: Root Node
    if (current_node == this->_root_node) {
      DEBUG_STDOUT("HANDLE_ROOT: Node is root");
      return DelCase::HANDLE_ROOT;
    }

    if (current_node_right == nullptr && current_node_left == nullptr) {
      DEBUG_STDOUT("HANDLE_ROOT: Left and right siblings do not exist");
      return DelCase::HANDLE_ROOT;
    }

    if (current_node_right != nullptr && current_node_left != nullptr) {
      DEBUG_STDOUT("COND 1: Both sibling nodes exist");
      // * Both siblings exist

      // ! Case 2: Right node has enough data and left too little
      if (current_node_right->is_above_min() &&
          current_node_left->is_below_or_at_min()) {
        if (parent_node == parent_node_right) {
          DEBUG_STDOUT("BALANCE_RIGHT_TO_CURRENT: Right node has enough data "
                       "and the same parent");
          return DelCase::BALANCE_RIGHT_TO_CURRENT;
        } else {
          DEBUG_STDOUT("MERGE_LEFT_AND_CURRENT: Left node would underflow and "
                       "can be merged");
          return DelCase::MERGE_LEFT_AND_CURRENT;
        }
      }

      // ! Case 3: Right node has too little data and left enough
      if (current_node_right->is_below_or_at_min() &&
          current_node_left->is_above_min()) {
        if (parent_node == parent_node_left) {
          DEBUG_STDOUT("BALANCE_LEFT_TO_CURRENT: Left node has enough data and "
                       "the same parent");
          return DelCase::BALANCE_LEFT_TO_CURRENT;
        } else {
          DEBUG_STDOUT("MERGE_RIGHT_AND_CURRENT: Right node would underflow "
                       "and can be merged");
          return DelCase::MERGE_RIGHT_AND_CURRENT;
        }
      }

      // ! Case 4: Both nodes have enough and the same parent
      // !   => Choose node that has more available
      if ((current_node_right->is_above_min() &&
           current_node_left->is_above_min()) &&
          (parent_node_left == parent_node_right)) {
        DEBUG_STDOUT(
            "Nodes have the same parent, choosing the larger to balance");
        if (current_node_right->size >= current_node_left->size) {
          DEBUG_STDOUT("BALANCE_RIGHT_TO_CURRENT: Right node is larger than "
                       "left node");
          return DelCase::BALANCE_RIGHT_TO_CURRENT;
        } else {
          DEBUG_STDOUT(
              "BALANCE_LEFT_TO_CURRENT: Left node is larger than right node");
          return DelCase::BALANCE_LEFT_TO_CURRENT;
        }
      }
    }

    if (current_node_right != nullptr && parent_node == parent_node_right) {
      DEBUG_STDOUT("COND 2: Right node exists and has same parent");

      // ! Case 5: Right node exists, has same parent, and can be merged
      if (current_node_right->is_below_or_at_min()) {
        DEBUG_STDOUT("Right node can be merged");
        return DelCase::MERGE_RIGHT_AND_CURRENT;
      }

      // ! Case 6: Right node exists, has same parent, and can be balanced
      DEBUG_STDOUT("Right node can be balanced");
      return DelCase::BALANCE_RIGHT_TO_CURRENT;
    }

    if (current_node_left != nullptr && parent_node == parent_node_left) {
      DEBUG_STDOUT("COND 3: Left node exists and has same parent");

      // ! Case 7: Left node exists, has same parent, and can be merged
      if (current_node_left->is_below_or_at_min()) {
        DEBUG_STDOUT("Left node can be merged");
        return DelCase::MERGE_LEFT_AND_CURRENT;
      }

      // ! Case 8: Left node exists, has same parent, and can be balanced
      DEBUG_STDOUT("Left Node can be balanced");
      return DelCase::BALANCE_LEFT_TO_CURRENT;
    }

    DEBUG_STDERR("COND 3: Neither sibling exists!");

    // ! Neither sibling exists; should NEVER happen!
    throw std::runtime_error("No siblings found - this should never happen!");
  }

  /*
  This method is genuinely painful.

  Proper deletion must happen recursively; the left and right nodes
  of the current iteration's node *and* its parent index node must
  be supplied for merging to work correctly.
  */
  DelState _erase(const key_type &key, Node *const &current_node,
                  Node *const &current_node_left = nullptr,
                  Node *const &current_node_right = nullptr,
                  IndexNode *const &parent_node = nullptr,
                  IndexNode *const &parent_node_left = nullptr,
                  IndexNode *const &parent_node_right = nullptr,
                  size_type *const &parent_position = nullptr) {

    DEBUG_STDOUT("-- Depth: " << (this->_root_node->layer - current_node->layer)
                              << " | Erasing " << key);

    if (current_node->is_leaf()) {
      LeafNode *l_node{dynamic_cast<LeafNode *>(current_node)};
      LeafNode *left_l_node{dynamic_cast<LeafNode *>(current_node_left)};
      LeafNode *right_l_node{dynamic_cast<LeafNode *>(current_node_right)};

      DEBUG_STDOUT("Encountered LeafNode:   " << l_node);
      DEBUG_STDOUT("Left Sibling LeafNode:  " << left_l_node);
      DEBUG_STDOUT("Right Sibling LeafNode: " << right_l_node);

      size_type delete_position{this->get_pos(l_node, key)};

      // 0: Return DelState::FAILURE if key was not found
      if (l_node->size <= delete_position ||
          !key_equal_to{}(l_node->values[delete_position], key)) {
        DEBUG_STDOUT("FAILURE: Key " << key << " was not found");
        return DelState::FAILURE;
      }

      DEBUG_STDOUT("Deleting [" << delete_position << "]: " << key);

      // 1: If key was found, move all keys after it one
      //    position forward and adjust leaf node size
      DEBUG_STDOUT("Shifting values forward");
      for (size_type i{delete_position + 1}; i < l_node->size; ++i) {
        DEBUG_STDOUT("values[" << i << "-" << 1 << "] = values[" << i << "] ("
                               << l_node->values[i] << ")");
        l_node->values[i - 1] = l_node->values[i];
      }
      l_node->size--;

      // 2: Handle underflow scenarios
      DelCase deletion_case{this->_determine_delcase(
          l_node, left_l_node, right_l_node, parent_node, parent_node_left,
          parent_node_right)};

      DEBUG_STDOUT("Handling LeafNode DelCase: " << deletion_case);

      switch (deletion_case) {
      case DelCase::HANDLE_ROOT:
        DEBUG_STDOUT("case DelCase::HANDLE_ROOT");
        break;

      case DelCase::NO_ACTION:
        DEBUG_STDOUT("case DelCase::NO_ACTION");
        break;

      case DelCase::BALANCE_RIGHT_TO_CURRENT:
        DEBUG_STDOUT("case DelCase::BALANCE_RIGHT_TO_CURRENT");

        this->balance_to_left(l_node, right_l_node, parent_node,
                              *parent_position);
        break;

      case DelCase::BALANCE_LEFT_TO_CURRENT:
        DEBUG_STDOUT("case DelCase::BALANCE_LEFT_TO_CURRENT");

        this->balance_to_right(left_l_node, l_node, parent_node,
                               *parent_position - 1);
        break;

      case DelCase::MERGE_LEFT_AND_CURRENT:
        DEBUG_STDOUT("case DelCase::MERGE_LEFT_AND_CURRENT");

        this->merge_nodes(left_l_node, l_node, parent_node,
                          *parent_position - 1);
        return DelState::MERGE;

      case DelCase::MERGE_RIGHT_AND_CURRENT:
        DEBUG_STDOUT("case DelCase::MERGE_RIGHT_AND_CURRENT");

        this->merge_nodes(l_node, right_l_node, parent_node, *parent_position);
        return DelState::MERGE;
      }

      DEBUG_STDOUT("SUCCESS: LeafNode cases should be handled");
      return DelState::SUCCESS;

    } else { // current_node != leaf
      IndexNode *i_node{dynamic_cast<IndexNode *>(current_node)};
      IndexNode *left_i_node{dynamic_cast<IndexNode *>(current_node_left)};
      IndexNode *right_i_node{dynamic_cast<IndexNode *>(current_node_right)};

      DEBUG_STDOUT("Encountered IndexNode:   " << i_node);
      DEBUG_STDOUT("Left Sibling IndexNode:  " << left_i_node);
      DEBUG_STDOUT("Right Sibling IndexNode: " << right_i_node);

      size_type delete_position{this->get_pos(i_node, key)};

      // ! If key at delete pos is equal to key, increment pos
      if (delete_position < i_node->size &&
          key_equal_to{}(key, i_node->at(delete_position))) {
        delete_position++;
      }

      DelState deletion_state{DelState::NON_INIT};

      while (delete_position < i_node->size + 1) {
        DEBUG_STDOUT("Deleting [" << delete_position << "]: " << key);

        Node *passed_current = i_node->children[delete_position];
        Node *passed_left = nullptr;
        Node *passed_right = nullptr;

        IndexNode *passed_parent = i_node;
        IndexNode *passed_parent_left = nullptr;
        IndexNode *passed_parent_right = nullptr;

        size_type *passed_position = &delete_position;

        // Lower and upper bounds of current node are checked here
        // and are never exceeded

        // Lower bound (left side / limit)
        if (delete_position == 0) {
          if (current_node_left == nullptr) {
            passed_left = nullptr;
          } else {
            // ! Last child node of left sibling node
            passed_left = left_i_node->children[left_i_node->size];
          }

        } else {
          // Pass along direct sibling to the left
          passed_left = i_node->children[delete_position - 1];
          passed_parent_left = i_node;
        }

        // Upper bound (right side / limit)
        if (delete_position == i_node->size) {
          if (current_node_right == nullptr) {
            passed_right = nullptr;
          } else {
            // ! First child node of right sibling node
            passed_right = right_i_node->children[0];
          }

        } else {
          // Pass along direct sibling to the right
          passed_right = i_node->children[delete_position + 1];
          passed_parent_right = i_node;
        }

        DEBUG_STDOUT("TRAVERSING DOWN");
        DEBUG_STDOUT("passed_current:      " << passed_current);
        DEBUG_STDOUT("passed_left:         " << passed_left);
        DEBUG_STDOUT("passed_right:        " << passed_right);
        DEBUG_STDOUT("passed_parent:       " << passed_parent);
        DEBUG_STDOUT("passed_parent_left:  " << passed_parent_left);
        DEBUG_STDOUT("passed_parent_right: " << passed_parent_right);
        DEBUG_STDOUT("passed_position:     " << *passed_position);

        deletion_state = this->_erase(
            key, passed_current, passed_left, passed_right, passed_parent,
            passed_parent_left, passed_parent_right, passed_position);

        if ((deletion_state == DelState::SUCCESS) ||
            (deletion_state == DelState::MERGE)) {
          DEBUG_STDOUT("Encountered DelState: " << deletion_state
                                                << " - breaking");
          break;
        }

        ++delete_position;
      }

      if (delete_position > i_node->size) {
        deletion_state = DelState::FAILURE;
      }

      DEBUG_STDOUT("DELETE: IndexNode after traversal: " << i_node);
      DEBUG_STDOUT("Handling IndexNode DelState: " << deletion_state);
      DEBUG_STDOUT("DELETE: Tree after traversal: ");
#ifdef DEBUG
      this->dump_detailed();
#endif

      // ! Deletion State Handling
      switch (deletion_state) {
      case DelState::FAILURE:
        DEBUG_STDOUT("case DelState::FAILURE");
        // So the compiler doesn't scream when -DDEBUG is set
        return deletion_state;

      case DelState::SUCCESS:
        DEBUG_STDOUT("case DelState::SUCCESS");
        return deletion_state;

      case DelState::MERGE:
        DEBUG_STDOUT("case DelState::MERGE");

        // ! Current or next node need to be shifted
        if (i_node->children[delete_position]->size != 0) {
          delete_position++;
        }

        DEBUG_STDOUT("\tdelete i_node->children[" << delete_position << "]");
        if (i_node->children[delete_position]->is_leaf()) {
          this->_l_node_cache->push(
              dynamic_cast<LeafNode *>(i_node->children[delete_position]));
          this->_leaf_nodes--;
        } else {
          this->_i_node_cache->push(
              dynamic_cast<IndexNode *>(i_node->children[delete_position]));
          this->_index_nodes--;
        }

        // ! delete_position must be decremented for keys
        // * merge_nodes() always merges right node to left,
        // * so delete_position should never be 0
        for (size_type i{delete_position - 1}; i < i_node->size - 1; ++i) {
          DEBUG_STDOUT("\ti_node->keys[" << i << "] = i_node->keys[" << i + 1
                                         << "]");
          i_node->keys[i] = i_node->keys[i + 1];
        }

        for (size_type i{delete_position}; i < i_node->size; ++i) {
          DEBUG_STDOUT("\ti_node->children[" << i << "] = i_node->children["
                                             << i + 1 << "] ("
                                             << i_node->children[i + 1] << ")");
          i_node->children[i] = i_node->children[i + 1];
        }

        DEBUG_STDOUT(
            "HANDLE MERGE 1: IndexNode after merge handle: " << i_node);

        i_node->size--;

        // this->repair_index_node(i_node);

        DEBUG_STDOUT("HANDLE MERGE 2: IndexNode after handling: ");
#ifdef DEBUG
        this->dump_detailed(std::cerr, i_node);
#endif
        break;

      case DelState::NON_INIT:
        DEBUG_STDOUT("case DelState::NON_INIT");

        throw std::runtime_error(
            "deletion_state was not initialized -- good job!");
      }

      // ! Deletion case handling
      DelCase deletion_case{this->_determine_delcase(
          i_node, left_i_node, right_i_node, parent_node, parent_node_left,
          parent_node_right)};

      DEBUG_STDOUT("Handling IndexNode DelCase: " << deletion_case);

      DelState new_deletion_state{DelState::NON_INIT};

      switch (deletion_case) {
      // ! If IndexNode is root and is empty, re-assign first leaf
      case DelCase::HANDLE_ROOT:
        DEBUG_STDOUT("case DelCase::HANDLE_ROOT");

        if (i_node->size == 0) {
          // * Root IndexNodes with size 0 are still expected to have
          // * exactly one child node left -- this node is then assigned
          // * as the new root
          this->_root_node = i_node->children[0];
          this->_i_node_cache->push(i_node);
          this->_index_nodes--;
          i_node = dynamic_cast<IndexNode *>(this->_root_node);

          DEBUG_STDOUT("HANDLE ROOT: Assigned "
                       << this->_root_node
                       << " as new root and deleted old node");
        }
        new_deletion_state = DelState::SUCCESS;
        break;

      case DelCase::BALANCE_RIGHT_TO_CURRENT:
        DEBUG_STDOUT("case DelCase::BALANCE_RIGHT_TO_CURRENT");

        this->balance_to_left(i_node, right_i_node, parent_node,
                              *parent_position);
        new_deletion_state = DelState::SUCCESS;
        break;

      case DelCase::BALANCE_LEFT_TO_CURRENT:
        DEBUG_STDOUT("case DelCase::BALANCE_LEFT_TO_CURRENT");

        this->balance_to_right(left_i_node, i_node, parent_node,
                               *parent_position - 1);
        new_deletion_state = DelState::SUCCESS;
        break;

      case DelCase::MERGE_LEFT_AND_CURRENT:
        DEBUG_STDOUT("case DelCase::MERGE_LEFT_AND_CURRENT");

        this->merge_nodes(left_i_node, i_node, parent_node,
                          *parent_position - 1);
        new_deletion_state = DelState::MERGE;
        break;

      case DelCase::MERGE_RIGHT_AND_CURRENT:
        DEBUG_STDOUT("case DelCase::MERGE_RIGHT_AND_CURRENT");

        this->merge_nodes(i_node, right_i_node, parent_node, *parent_position);
        new_deletion_state = DelState::MERGE;
        break;

      case DelCase::NO_ACTION:
        DEBUG_STDOUT("case DelCase::NO_ACTION");
        new_deletion_state = DelState::SUCCESS;
        break;
      }

      DEBUG_STDOUT("Final IndexNode: " << i_node);
      DEBUG_STDOUT("SUCCESS: IndexNode cases should be handled");
      return new_deletion_state;
    }
  }

  /*
  TODO: Docstring
  */
  void merge_nodes(Node *const &left_node, Node *const &right_node,
                   IndexNode *const &parent_node = nullptr,
                   const size_type &parent_key_pos = size_type()) {
    if (left_node->is_leaf() != right_node->is_leaf()) {
      throw std::runtime_error("Cannot merge nodes of mismatching type");
    }

    if (left_node->size + right_node->size > Properties::max_store) {
      throw std::runtime_error("Node merge would cause overflow");
    }

    DEBUG_STDOUT(
        "MERGE: " << (left_node->is_leaf() ? "LeafNode" : "IndexNode"));
    DEBUG_STDOUT("\tPosition: " << parent_key_pos);
    DEBUG_STDOUT("\tBefore:");
    DEBUG_STDOUT("\tLeft:  " << left_node);
    DEBUG_STDOUT("\tRight: " << right_node);
    DEBUG_STDOUT("\tParent: " << parent_node);

    if (left_node->is_leaf()) {

      LeafNode *left_l_node{dynamic_cast<LeafNode *>(left_node)};
      LeafNode *right_l_node{dynamic_cast<LeafNode *>(right_node)};

      // 1: Copy values from right to left
      for (size_type i{0}; i < right_l_node->size; ++i) {
        left_l_node->values[left_l_node->size + i] = right_l_node->values[i];
      }

      // 3: Adjust links of left leaf
      left_l_node->right_leaf = right_l_node->right_leaf;
      if (left_l_node->right_leaf != nullptr) {
        left_l_node->right_leaf->left_leaf = left_l_node;
      } else {
        this->_rightmost_leaf = left_l_node;
      }

      // 4: Update parent node
      // ! Special case: left_l_node->size != 0 for N=1
      // ! Parent Update
      if (parent_key_pos > 0 && left_l_node->size != 0) {
        DEBUG_STDOUT("MERGE: Updating parent node at position "
                     << parent_key_pos - 1 << " with first key of left node: "
                     << right_l_node->at(0));
        parent_node->keys[parent_key_pos - 1] = left_l_node->at(0);
      } else {
        DEBUG_STDOUT("MERGE: Not updating parent node at position "
                     << parent_key_pos << " - 1");
      }

    } else {
      // ! Keep in mind that IndexNode merges imply that there's
      // ! a parent node, so we're not checking for that here

      IndexNode *left_i_node{dynamic_cast<IndexNode *>(left_node)};
      IndexNode *right_i_node{dynamic_cast<IndexNode *>(right_node)};

      // 1: Add parent key
      // The parent's key is necessary because n keys,
      // but n+1 children are moved
      DEBUG_STDOUT("\tParent key: " << parent_node->keys[parent_key_pos]);
      left_i_node->keys[left_i_node->size] = parent_node->keys[parent_key_pos];
      left_i_node->size++;

      // 2: Copy keys and children from right to left
      for (size_type i{0}; i < right_i_node->size; ++i) {
        left_i_node->keys[left_i_node->size + i] = right_i_node->keys[i];
        left_i_node->children[left_i_node->size + i] =
            right_i_node->children[i];
      }

      // 3: Append last child
      left_i_node->children[left_i_node->size + right_i_node->size] =
          right_i_node->children[right_i_node->size];
    }

    // Adjust size
    left_node->size += right_node->size;
    right_node->size = 0;

    DEBUG_STDOUT("\tAfter: ");
    DEBUG_STDOUT("\tLeft:  " << left_node);
    DEBUG_STDOUT("\tRight: " << right_node);
  }

  /*
  ! DEPRECATED: turns out I'm just stupid

  For some reason that I cannot possibly comprehend,
  there seems to be an issue with index nodes allowing
  duplicate inserts. My gut feeling tells me that it
  happens in very specific cases during node merges.

  This is a very poor attempt at a fix.
  */
  size_type repair_index_node(IndexNode *&i_node) {
    DEBUG_STDOUT("REPAIR: Checking for misaligned IndexNode keys");
    DEBUG_STDOUT("index node: " << i_node);
    if (i_node->size == 0) {
      return 0;
    }

    size_type count{0};

    auto child_key_getter = [](Node *n) -> key_type {
      if (n->is_leaf()) {
        auto l{dynamic_cast<LeafNode *>(n)};
        return l->values[0];
      } else {
        auto i{dynamic_cast<IndexNode *>(n)};
        return i->keys[0];
      }
    };

    for (size_type i{0}; i < i_node->size; ++i) {
      key_type child_key = child_key_getter(i_node->children[i + 1]);
      DEBUG_STDOUT("\tchild node: " << i_node->children[i + 1]);
      DEBUG_STDOUT("\tchild_key = " << child_key);

      if (key_compare{}(child_key, i_node->keys[i])) {
        DEBUG_STDOUT("\tFIXUP INDEX KEY: i_node->keys["
                     << i << "] (" << i_node->keys[i] << ") = " << child_key);
        i_node->keys[i] = child_key;
        count++;
      }
    }

    DEBUG_STDOUT("REPAIR END: Fixed " << count << " keys in total.");
    return count;
  }

  /*
  TODO: Docstring
  */
  void balance_to_left(Node *const &left_node, Node *const &right_node,
                       IndexNode *const &parent_node = nullptr,
                       const size_type &parent_update_pos = size_type()) {
    // 0: Ensure node types are equal
    if (left_node->is_leaf() != right_node->is_leaf()) {
      throw std::runtime_error(
          "Cannot balance nodes of mismatching type to the left");
    }

    // 1: Determine number of values to move => half of total size
    size_type move_count{(right_node->size - left_node->size) / 2};

    if (left_node->is_leaf()) {

      LeafNode *left_l_node{dynamic_cast<LeafNode *>(left_node)};
      LeafNode *right_l_node{dynamic_cast<LeafNode *>(right_node)};

      DEBUG_STDOUT("BALANCE_TO_LEFT: LeafNodes");
      DEBUG_STDOUT("\tPosition: " << parent_update_pos);
      DEBUG_STDOUT("\tBefore:");
      DEBUG_STDOUT("\tLeft:  " << left_l_node);
      DEBUG_STDOUT("\tRight: " << right_l_node);

      // 2: Append values from right node to left node
      for (size_type i{0}; i < move_count; ++i) {
        left_l_node->values[left_l_node->size + i] = right_l_node->values[i];
      }

      // 3: Shift remaining values of right size forward
      for (size_type i{0}; i < right_l_node->size - move_count; ++i) {
        right_l_node->values[i] = right_l_node->values[i + move_count];
      }

      // 4: Adjust node sizes
      left_l_node->size += move_count;
      right_l_node->size -= move_count;

      // 5: Set key of parent node to first element of right node
      // ! Parent Update
      parent_node->keys[parent_update_pos] = right_l_node->values[0];
      DEBUG_STDOUT("\t - parent_node->keys["
                   << parent_update_pos << "] = " << right_l_node->values[0]);

      DEBUG_STDOUT("\tAfter: ");
      DEBUG_STDOUT("\tLeft:  " << left_l_node);
      DEBUG_STDOUT("\tRight: " << right_l_node);

    } else {

      // 2: Ensure parent node exists
      if (parent_node == nullptr) {
        throw std::runtime_error(
            "No parent node found while balancing index nodes to left");
      }

      IndexNode *left_i_node{dynamic_cast<IndexNode *>(left_node)};
      IndexNode *right_i_node{dynamic_cast<IndexNode *>(right_node)};

      DEBUG_STDOUT("BALANCE_TO_LEFT: IndexNodes");
      DEBUG_STDOUT("\tLeft:  " << left_i_node);
      DEBUG_STDOUT("\tRight: " << right_i_node);
      DEBUG_STDOUT("\tParent Position: " << parent_update_pos);

      // 3: Append key of parent to left node => "inbetween" key
      // because keys until move_count - 1 are shifted, but
      // child nodes until move_count
      DEBUG_STDOUT("\nBALANCE_TO_LEFT: Append key from parent");
      DEBUG_STDOUT("\tleft_i_node->keys["
                   << left_i_node->size << "] = parent_node->keys["
                   << parent_update_pos << "] ("
                   << parent_node->keys[parent_update_pos] << ")");
      DEBUG_STDOUT("\tleft_i_node->size + 1 = " << left_i_node->size + 1);

      left_i_node->keys[left_i_node->size] =
          parent_node->keys[parent_update_pos];
      left_i_node->size++;

      // 4: Append keys from right to left node
      DEBUG_STDOUT("\nBALANCE_TO_LEFT: Appending keys from right to left node");
      for (size_type i{0}, limit{move_count - 1}; i < limit; ++i) {
        DEBUG_STDOUT("\tleft_i_node->keys[" << left_i_node->size + i
                                            << "] = right_i_node->keys[" << i
                                            << "]");

        left_i_node->keys[left_i_node->size + i] = right_i_node->keys[i];
      }

      // 5: Append children from right to left node
      DEBUG_STDOUT(
          "\nBALANCE_TO_LEFT: Appending children from right to left node");
      for (size_type i{0}; i < move_count; ++i) {
        DEBUG_STDOUT("\tleft_i_node->children[" << left_i_node->size + i
                                                << "] = right_i_node->children["
                                                << i << "]");

        left_i_node->children[left_i_node->size + i] =
            right_i_node->children[i];
      }

      /*
      for (size_type i{0}; i < move_count; ++i) {
        DEBUG_STDOUT("\tleft_i_node->keys[" << left_i_node->size + i
                                            << "] = right_i_node->keys[" << i
                                            << "]");
        DEBUG_STDOUT("\tleft_i_node->children[" << left_i_node->size + i
                                                << "] =
      right_i_node->children["
                                                << i << "]");

        left_i_node->keys[left_i_node->size + i] = right_i_node->keys[i];
        left_i_node->children[left_i_node->size + i] =
            right_i_node->children[i];
      }

      // 5: Append last child
      DEBUG_STDOUT("\nBALANCE_TO_LEFT: Appending last child node");
      DEBUG_STDOUT("\tleft_i_node->children[" << left_i_node->size +
      move_count
                                              << "] = right_i_node->children["
                                              << move_count << "]");

      left_i_node->children[left_i_node->size + move_count] =
          right_i_node->children[move_count];
      */

      // 6: Update size of left node
      DEBUG_STDOUT("\nBALANCE_TO_LEFT: Updating size of left node");
      DEBUG_STDOUT("left_i_node->size + move_count - 1 = "
                   << left_i_node->size + move_count - 1);

      left_i_node->size += move_count - 1;

      // 7: Update parent index node's key with the first
      // key after those that were appended -- that's "the key in the middle"
      // ! this key will be overwritten in the next step
      // ! Parent Update
      parent_node->keys[parent_update_pos] = right_i_node->keys[move_count - 1];
      DEBUG_STDOUT("\nBALANCE_TO_LEFT: Updating parent index node's key");
      DEBUG_STDOUT("\t - parent_node->keys["
                   << parent_update_pos
                   << "] = " << right_i_node->keys[move_count - 1]);

      /*
      This is hot garbage but could maybe used for optimizations
      // 8: Shift remaining keys and child nodes of right node forward
      for (size_type i{move_count}; i < right_i_node->size; ++i) {
        right_i_node->keys[i - move_count] = right_i_node->keys[i];
        right_i_node->children[i - move_count] = right_i_node->children[i];
      }

      // 9: Shift remaining child node forward
      right_i_node->children[move_count + 1] =
          right_i_node->children[right_i_node->size];
      */

      // 8: Move keys of right node forward
      DEBUG_STDOUT("\nBALANCE_TO_LEFT: Move keys of right node forward");

      for (size_type i{0}, limit{right_i_node->size - move_count}; i < limit;
           ++i) {
        DEBUG_STDOUT("\tright_i_node->keys[" << i << "] = right_i_node->keys["
                                             << i + move_count << "]");

        right_i_node->keys[i] = right_i_node->keys[i + move_count];
      }

      // 9: Move child nodes of right node forward
      // * See the + 1 there? Took me so long to figure that out
      DEBUG_STDOUT("\nBALANCE_TO_LEFT: Move children of right node forward");

      for (size_type i{0}, limit{right_i_node->size - move_count + 1};
           i < limit; ++i) {
        DEBUG_STDOUT("\tright_i_node->children["
                     << i << "] = right_i_node->children[" << i + move_count
                     << "]");

        right_i_node->children[i] = right_i_node->children[i + move_count];
      }

      // 10: Finally, decrease size of right node
      DEBUG_STDOUT("\nBALANCE_TO_LEFT: Decrease size of right node");
      DEBUG_STDOUT("right_i_node->size - move_count = " << right_i_node->size -
                                                               move_count);

      right_i_node->size -= move_count;

      DEBUG_STDOUT("\nBALANCE_TO_LEFT: RESULT");
      DEBUG_STDOUT("\tLeft:  " << left_i_node);
      DEBUG_STDOUT("\tRight: " << right_i_node);
    }
  }

  /*
  TODO: Docstring
  */
  void balance_to_right(Node *const &left_node, Node *const &right_node,
                        IndexNode *const &parent_node = nullptr,
                        const size_type &parent_update_pos = size_type()) {
    // 0: Ensure node types are equal
    if (left_node->is_leaf() != right_node->is_leaf()) {
      throw std::runtime_error(
          "Cannot balance nodes of mismatching type to the right");
    }

    // 1: Determine number of values to move => half of total size
    size_type move_count{(left_node->size - right_node->size) / 2};

    if (left_node->is_leaf()) {

      LeafNode *left_l_node{dynamic_cast<LeafNode *>(left_node)};
      LeafNode *right_l_node{dynamic_cast<LeafNode *>(right_node)};

      DEBUG_STDOUT("BALANCE_TO_RIGHT: LeafNodes");
      DEBUG_STDOUT("\tPosition: " << parent_update_pos);
      DEBUG_STDOUT("\tBefore:");
      DEBUG_STDOUT("\tLeft:  " << left_l_node);
      DEBUG_STDOUT("\tRight: " << right_l_node);

      // 2: Shift the values of the right node further back
      for (size_type i{right_l_node->size}; i-- > 0;) {
        right_l_node->values[i + move_count] = right_l_node->values[i];
      }

      // 3: Prepend values from left node to right node
      for (size_type l{left_l_node->size - move_count}, r{0};
           l < left_l_node->size; ++l) {
        right_l_node->values[r++] = left_l_node->values[l];
      }

      // 4: Adjust node sizes
      right_l_node->size += move_count;
      left_l_node->size -= move_count;

      // 5: Set key of parent node to first element of right node
      // ! Parent Update
      parent_node->keys[parent_update_pos] = right_l_node->values[0];
      DEBUG_STDOUT("\t - parent_node->keys["
                   << parent_update_pos << "] = " << right_l_node->values[0]);

      DEBUG_STDOUT("\tAfter: ");
      DEBUG_STDOUT("\tLeft:  " << left_l_node);
      DEBUG_STDOUT("\tRight: " << right_l_node);

    } else {
      // ! ! ! Note: This gave me an aneurysm; DO NOT TOUCH

      // 2: Ensure parent node exists
      if (parent_node == nullptr) {
        throw std::runtime_error(
            "No parent node found while balancing index nodes to left");
      }

      IndexNode *left_i_node{dynamic_cast<IndexNode *>(left_node)};
      IndexNode *right_i_node{dynamic_cast<IndexNode *>(right_node)};

      DEBUG_STDOUT("BALANCE_TO_RIGHT: IndexNodes");
      DEBUG_STDOUT("\tPosition: " << parent_update_pos);
      DEBUG_STDOUT("\tBefore:");
      DEBUG_STDOUT("\tLeft:  " << left_i_node);
      DEBUG_STDOUT("\tRight: " << right_i_node);

      // 3: Shift keys and children to the right
      DEBUG_STDOUT("\nBALANCE_TO_RIGHT: Shift keys and children to the right");
      for (size_type i{right_i_node->size}; i-- > 0;) {
        DEBUG_STDOUT("\tright_i_node->keys[" << i << " + " << move_count
                                             << "] = right_i_node->keys[" << i
                                             << "]");
        DEBUG_STDOUT("\tright_i_node->children["
                     << i << " + 1 + " << move_count
                     << "] = right_i_node->children[" << i << " + 1]");

        right_i_node->keys[i + move_count] = right_i_node->keys[i];
        right_i_node->children[i + 1 + move_count] =
            right_i_node->children[i + 1];
      }

      // 4: Shift remaining child at the beginning
      DEBUG_STDOUT("\nBALANCE_TO_RIGHT: Shift remaining child at beginning");
      DEBUG_STDOUT("\tright_i_node->children["
                   << move_count << "] = right_i_node->children[0]");

      right_i_node->children[move_count] = right_i_node->children[0];

      // 5: Insert key of parent node
      DEBUG_STDOUT("\nBALANCE_TO_RIGHT: Insert key of parent node");
      DEBUG_STDOUT("\tright_i_node->keys["
                   << move_count << " - 1] = parent_node->keys["
                   << parent_update_pos << "] ("
                   << parent_node->keys[parent_update_pos] << ")");

      right_i_node->keys[move_count - 1] = parent_node->keys[parent_update_pos];

      // ! 6: ... and update the key between the two nodes in the parent node
      // ! This here again is the "middle key" and was determined through
      // ! trial and error.
      // ! Note to self: Don't touch it please please please please
      // ! Parent Update
      DEBUG_STDOUT("\nBALANCE_TO_RIGHT: Insert middle key in parent parent");
      DEBUG_STDOUT("\tparent_node->keys["
                   << parent_update_pos << "] = "
                   << left_i_node->keys[left_i_node->size - move_count]);

      parent_node->keys[parent_update_pos] =
          left_i_node->keys[left_i_node->size - move_count];

      // 7: Update size of right node
      DEBUG_STDOUT("\nBALANCE_TO_RIGHT: Update size of right node; old: "
                   << right_i_node->size
                   << "; new: " << right_i_node->size + move_count);

      right_i_node->size += move_count;

      // 8: Prepend keys from left to right
      // * move_count - 1 because we already prepended the parent node's key
      // * with the middle key before
      DEBUG_STDOUT("\nBALANCE_TO_RIGHT: Prepending keys from left to right");

      for (size_type r{move_count - 1}, l{left_i_node->size - 1}; r-- > 0;
           --l) {
        DEBUG_STDOUT("\tright_i_node->[" << r << "] = left_i_node->[" << l
                                         << "]");

        right_i_node->keys[r] = left_i_node->keys[l];
      }

      // 9: Prepend children from left to right
      DEBUG_STDOUT(
          "\nBALANCE_TO_RIGHT: Prepending children from left to right");
      for (size_type r{move_count}, l{left_i_node->size}; r-- > 0; --l) {
        DEBUG_STDOUT("\tright_i_node->[" << r << "] = left_i_node->[" << l
                                         << "]");

        right_i_node->children[r] = left_i_node->children[l];
      }

      /*
      // Hot garbage that doesn't handle move_count == 1
      * Even though this is shit, this could be used for optimization
      purposes? for (size_t l{left_i_node->size - move_count}, r{0}; l <
      left_i_node->size;) { DEBUG_STDOUT("\tright_i_node->keys[" << r << "] =
      left_i_node->keys["
                                             << l << "] ("
                                             << left_i_node->keys[l] << ")");
        DEBUG_STDOUT("\tright_i_node->children["
                     << r << "] = left_i_node->children[" << l + 1 << "] ("
                     << left_i_node->children[l + 1] << ")");

        right_i_node->keys[r] = left_i_node->keys[l];
        right_i_node->children[r++] = left_i_node->children[++l];
      }
      */

      // 10: Update size of left node, finally
      DEBUG_STDOUT("\nBALANCE_TO_RIGHT: Update size of left node; old: "
                   << left_i_node->size
                   << "; new: " << left_i_node->size - move_count);

      left_i_node->size -= move_count;

      DEBUG_STDOUT("\n\tAfter: ");
      DEBUG_STDOUT("\tLeft:  " << left_i_node);
      DEBUG_STDOUT("\tRight: " << right_i_node);
    }
  }

  // * * * Iterator Implementation * * *
public:
  class ConstIterator {
  public:
    using key_type = bplustree::key_type; // Do I need this?
    using value_type = bplustree::value_type;
    using reference = bplustree::const_reference;
    using pointer = bplustree::const_pointer;
    using difference_type = bplustree::difference_type;
    using iterator_category = std::forward_iterator_tag;

  private:
    friend class bplustree<key_type, N, MaxCacheSize, Properties>;

    LeafNode *l_node;
    size_type current_pos;

    /*
    Does what it says on the tin. Advances the iterator.
    Added for DRY's sake.
    */
    void advance() {
      if (this->current_pos < this->l_node->size) {
        this->current_pos++;

        if (this->current_pos == this->l_node->size &&
            this->l_node->right_leaf != nullptr) {
          this->current_pos = 0;
          this->l_node = this->l_node->right_leaf;
        }
      }
    }

  public:
    // * * Constructors * *
    ConstIterator() : l_node{nullptr}, current_pos{0} {}

    ConstIterator(LeafNode *l_node) : l_node{l_node}, current_pos{0} {}
    ConstIterator(const LeafNode *l_node) : l_node{l_node}, current_pos{0} {}

    ConstIterator(LeafNode *l_node, const size_type &current_pos)
        : l_node{l_node}, current_pos{current_pos} {}
    ConstIterator(const LeafNode *l_node, const size_type &current_pos)
        : l_node{l_node}, current_pos{current_pos} {}

    // * * Operators * *
    reference operator*() const {
      return this->l_node->values[this->current_pos];
    }

    pointer operator->() const {
      return &this->l_node->values[this->current_pos];
    }

    bool operator==(const const_iterator &other) const {
      return (this->l_node == other.l_node &&
              this->current_pos == other.current_pos);
    }

    bool operator!=(const const_iterator &other) const {
      return !(*this == other);
    }

    const_iterator &operator++() { // prefix
      this->advance();
      return *this;
    }

    const_iterator operator++(int) { // postfix
      const_iterator initial = *this;
      this->advance();
      return initial;
    }
  };
};

template <typename Key, size_t N>
bool operator==(const bplustree<Key, N> &left, const bplustree<Key, N> &right) {
  if (left.size() != right.size()) {
    return false;
  }

  auto l_begin{left.begin()}, l_end{left.end()};
  auto r_begin{right.begin()}, r_end{right.end()};

  while (l_begin != l_end && r_begin != r_end) {
    if (!std::equal_to<Key>{}(*l_begin++, *r_begin++)) {
      return false;
    }
  }

  return true;
}

template <typename Key, size_t N>
bool operator!=(const zra::bplustree<Key, N> &left,
                const zra::bplustree<Key, N> &right) {
  return !(left == right);
}

template <typename Key, size_t N>
void swap(zra::bplustree<Key, N> &left, zra::bplustree<Key, N> &right) {
  left.swap(right);
}

} // namespace zra

#endif // BPLUSTREE_H