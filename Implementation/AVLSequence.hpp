/**
 * An implementation of the intrusive AVL tree to store a sequence of elements.
 * The inorder traversal of the tree coincides with the sequence.
 */

#pragma once

#include <iterator>
#include <string>
#include <iostream>
#include <algorithm>
#include <utility>
#include <functional>
#include <type_traits>
#include <boost/optional.hpp>

namespace Snu {
namespace Cnrc {
namespace AVLSequence {

class AVLSequenceModule {
	
public:
	
	template<class T, class SizeType>
	class AVLSequenceNodeMixin;
		
	/**
	 * This trait defines an interface for the type `T` to be stored in an intrusive AVL tree.
	 * If `T` extends AVLSequenceNodeMixin<T>, the default implementation of the trait will just work.
	 * If not, a specialization of this template for `T` must be provided.
	 */
	template<class T>
	class AVLSequenceNodeTrait {
	
	private:
		
		using NodeMixin = AVLSequenceNodeMixin<T, typename T::Size>;
		
		static_assert(
			std::is_base_of<NodeMixin, T>::value, 
			"The default implementation of AVLSequenceNodeTrait<T> requires T extending AVLSequenceNodeMixin. "
			"Or you can specialize AVLSequenceNodeTrait for T"
		);
		
	public:
		
		using Node = T;
		using Size = typename Node::Size;
		using Height = typename Node::Height;
	
	public:
		
		static T* parent(const T* const n) {
			return static_cast<const NodeMixin*>(n)->parent;
		}
		
		static void parent(T* const n, T* const p) {
			static_cast<NodeMixin*>(n)->parent = p;
		}
		
		static T* leftChild(const T* const n) {
			return static_cast<const NodeMixin*>(n)->leftChild;
		}
		
		static void leftChild(T* const n, T* const l) {
			static_cast<NodeMixin*>(n)->leftChild = l;
		}
	
		static T* rightChild(const T* n) {
			return static_cast<const NodeMixin*>(n)->rightChild;
		}
		
		static void rightChild(T* const n, T* const r) {
			static_cast<NodeMixin*>(n)->rightChild = r;
		}
		
		static Height height(const T* const n) {
			return static_cast<const NodeMixin*>(n)->height;
		}
		
		static void height(T* const n, Height h) {
			static_cast<NodeMixin*>(n)->height = h;
		}
		
		static Size volume(const T* const n) {
			return static_cast<const NodeMixin*>(n)->volume;
		}
		
		static void volume(T* const n, Size v) {
			static_cast<NodeMixin*>(n)->volume = v;
		}
		
		static void augment(T* const) {
		}
		
	};
	
	/**
	 * This class works as a mixin by CRTP to store a sequence of type `T` objects in 
	 * an AVL tree with the help of `AVLSequenceNodeTrait<T>`.
	 */
	template<class T, class SizeType>
	class AVLSequenceNodeMixin {
		
	public:
		
		using Size = SizeType;
		using Height = unsigned int;
	
	public:
		
		AVLSequenceNodeMixin()
		: parent(nullptr), leftChild(nullptr), rightChild(nullptr)
		, height(1), volume(1) {
		}
		
		AVLSequenceNodeMixin(const AVLSequenceNodeMixin&) = delete;
		
		AVLSequenceNodeMixin(AVLSequenceNodeMixin&& other) {
			*this = std::move(other);
		}
		
		AVLSequenceNodeMixin& operator=(const AVLSequenceNodeMixin&) = delete;
		
		AVLSequenceNodeMixin& operator=(AVLSequenceNodeMixin&& other) {
			parent = other.parent;
			leftChild = other.leftChild;
			rightChild = other.rightChild;
			height = other.height;
			volume = other.volume;
			if(parent != nullptr) {
				if(parent->leftChild == &other) {
					parent->leftChild = static_cast<T*>(this);
				} else {
					parent->rightChild = static_cast<T*>(this);
				}
			}
			if(leftChild != nullptr) {
				leftChild->parent = static_cast<T*>(this);
			}
			if(rightChild != nullptr) {
				rightChild->parent = static_cast<T*>(this);
			}
			other.parent = nullptr;
			other.leftChild = nullptr;
			other.rightChild = nullptr;
			other.height = 1;
			other.volume = 1;
			return *this;
		}
	
	private:
		
		T* parent;
		T* leftChild;
		T* rightChild;
		Height height;
		Size volume;

		friend class AVLSequenceNodeTrait<T>;
		
	};
	
	/**
	 * This provides operations to store a sequence of instances of `Node` 
	 * in an intrusive AVL tree with the aid of `AVLSequenceNodeTrait<Node>`.
	 * `Node` must implement the interface defined by the template class `AVLSequenceNodeTrait<Node>`.
	 * Or you can provide an equivalent trait implementation by specifying the second template paramater.
	 */
	template<class Node, class NodeTrait = AVLSequenceNodeTrait<Node>>
	class AVLSequenceAlgorithm {
		
	public:
		
		using Size = typename NodeTrait::Size;
		using Height = typename NodeTrait::Height;

	private:

		template<class T>
		using optional = boost::optional<std::reference_wrapper<T>>;
		
	private:
		
		using T = NodeTrait;
		
	private:
	
		template<class CNode>
		class IteratorImpl : public std::iterator<std::bidirectional_iterator_tag, CNode> {
			
		private:
			
			using This = IteratorImpl<CNode>;
		
		public:
			
			IteratorImpl() 
			: p(nullptr), q(nullptr) {
			}
			
			IteratorImpl(CNode* p, CNode* q)
			: p(p), q(q) {
			}
			
		public:
			
			This& operator++() {
				if(! p) {
					p = q;
					q = nullptr;
				} else if(auto r = next(p)) {
					p = &(r->get());
				} else {
					q = p;
					p = nullptr;
				}
				return *this;
			}
			
			This operator++(int) {
				This tmp(*this);
				operator++();
				return tmp;
			}
			
			This& operator--() {
				if(! p) {
					p = q;
					q = nullptr;
				} else if(auto r = previous(p)) {
					p = &(r->get());
				} else {
					q = p;
					p = nullptr;
				}
				return *this;
			}
			
			This operator--(int) {
				This tmp(*this);
				operator--();
				return tmp;
			}
			
			bool operator==(const This& other) const {
				return p == other.p && q == other.q;
			}
			
			bool operator!=(const This& other) const {
				return p != other.p || q != other.q;
			}
			
			CNode& operator*() const {
				return *p;
			}
			
			CNode* operator->() const {
				return p;
			}
		
		private:
			
			CNode* p;
			CNode* q;
		};
		
	public:
		
		class ContainerView {
			
		public:
			
			using iterator = IteratorImpl<Node>;
			using const_iterator = IteratorImpl<const Node>;
			using reverse_iterator = std::reverse_iterator<iterator>;
			using const_reverse_iterator = std::reverse_iterator<const_iterator>;
			
		public:
			
			explicit ContainerView(Node& root)
			: root(&root) {
			}
			
		public:
			
			iterator begin() {
				return iterator(findHead(root), nullptr);
			}
			
			iterator end() {
				return iterator(nullptr, findTail(root));
			}
			
			const_iterator begin() const {
				return const_iterator(findHead(root), nullptr);
			}
			
			const_iterator end() const {
				return const_iterator(nullptr, findTail(root));
			}
			
			reverse_iterator rbegin() {
				return std::reverse_iterator<iterator>(end());
			}
			
			reverse_iterator rend() {
				return std::reverse_iterator<iterator>(begin());
			}
			
			const_reverse_iterator rbegin() const {
				return std::reverse_iterator<const_iterator>(end());
			}
			
			const_reverse_iterator rend() const {
				return std::reverse_iterator<const_iterator>(begin());
			}
			
			const_iterator cbegin() const {
				return const_iterator(findHead(root), nullptr);
			}
			
			const_iterator cend() const {
				return const_iterator(nullptr, findTail(root));
			}
			
			const_reverse_iterator crbegin() const {
				return std::reverse_iterator<const_iterator>(end());
			}
			
			const_reverse_iterator crend() const {
				return std::reverse_iterator<const_iterator>(begin());
			}
			
		private:
			
			Node* const root;
		};
			
		class ConstContainerView {
			
		public:
			
			using const_iterator = IteratorImpl<const Node>;
			using const_reverse_iterator = std::reverse_iterator<const_iterator>;
			
		public:
			
			explicit ConstContainerView(const Node& root)
			: root(&root) {
			}
			
		public:
			
			const_iterator begin() const {
				return const_iterator(findHead(root), nullptr);
			}
			
			const_iterator end() const {
				return const_iterator(nullptr, findTail(root));
			}
			
			const_reverse_iterator rbegin() const {
				return std::reverse_iterator<const_iterator>(end());
			}
			
			const_reverse_iterator rend() const {
				return std::reverse_iterator<const_iterator>(begin());
			}
			
			const_iterator cbegin() const {
				return const_iterator(findHead(root), nullptr);
			}
			
			const_iterator cend() const {
				return const_iterator(nullptr, findTail(root));
			}
			
			const_reverse_iterator crbegin() const {
				return std::reverse_iterator<const_iterator>(end());
			}
			
			const_reverse_iterator crend() const {
				return std::reverse_iterator<const_iterator>(begin());
			}
			
		private:
			
			const Node* const root;
		};
		
	public:
		
		/**
		 * Returns ContainerView of the sequence whose root is `p`.
		 * ContainerView has an interface for STL compatible iterators.
		 */
		static ContainerView containerView(Node& p) {
			return ContainerView(p);
		}
		
		/**
		 * Returns ContainerView of the sequence whose root is `p`.
		 * ContainerView has an interface for STL compatible iterators.
		 */
		static ConstContainerView containerView(const Node& p) {
			return ConstContainerView(p);
		}
		
		/**
		 * Inserts node `n` before `p`.
		 */
		static void insertNodeBefore(Node& p, Node& n) {
			insertNodeBefore(&p, &n);
		}
		
		/**
		 * Inserts node `n` after `p`.
		 */
		static void insertNodeAfter(Node& p, Node& n) {
			insertNodeAfter(&p, &n);
		}
		
		/**
		 * Isolate `p` from its tree.
		 * `p` becomes a tree of its own.
		 */
		static void removeNode(Node& p) {
			removeNode(&p);
		}
		
		/**
		 * Joins the sequence whose rightmost element is `p` with 
		 * the sequence whose leftmost element is `q`.
		 */
		static void join(Node& p, Node& q) {
			join(&p, &q);
		}
		
		/**
		 * Splits `p` and its predecessor.
		 */
		static void splitBefore(Node& p) {
			splitBefore(&p);
		}
		
		/**
		 * Splits `p` and its successor.
		 */
		static void splitAfter(Node& p) {
			splitAfter(&p);
		}
		
		/**
		 * Returns the size of the subtree rooted at `p`.
		 */
		static Size size(const Node& p) {
			return T::volume(&p);
		}
		
		/**
		 * Returns the root of the tree containing `p`.
		 */
		static Node& findRoot(Node& p) {
			return *findRoot(&p); 
		}
		
		/**
		 * Returns the root of the tree containing `p`.
		 */
		static const Node& findRoot(const Node& p) {
			return *findRoot(&p); 
		}
		
		/**
		 * Finds the leftmost node of the subtree rooted at `p`.
		 * Call this with the root to get the first element of the whole sequence.
		 */
		static Node& findHead(Node& p) {
			return *findHead(&p);
		}
		
		/**
		 * Finds the leftmost node of the subtree rooted at `p`.
		 * Call this with the root to get the first element of the whole sequence.
		 */
		static const Node& findHead(const Node& p) {
			return *findHead(&p);
		}
		
		/**
		 * Finds the rightmost node of the subtree rooted at `p`.
		 * Call this with the root to get the last element of the whole sequence.
		 */
		static Node& findTail(Node& p) {
			return *findTail(&p);
		}
		
		/**
		 * Finds the rightmost node of the subtree rooted at `p`.
		 * Call this with the root to get the last element of the whole sequence.
		 */
		static const Node& findTail(const Node& p) {
			return *findTail(&p);
		}
		
		/**
		 * Returns the node next to `p` if `p` is not the last element.
		 */
		static optional<Node> next(Node& p) {
			return next(&p);
		}
		
		/**
		 * Returns the node next to `p` if `p` is not the last element.
		 */
		static optional<const Node> next(const Node& p) {
			return next(&p);
		}
		
		/**
		 * Returns the node previous to `p` if `p` is not the first element.
		 */
		static optional<Node> previous(Node& p) {
			return previous(&p);
		}
		
		/**
		 * Returns the node previous to `p` if `p` is not the first element.
		 */
		static optional<const Node> previous(const Node& p) {
			return previous(&p);
		}
	
	private:
		
		static void transplant(Node* p, Node* q) {
			if(T::parent(p)) {
				if(p == T::leftChild(T::parent(p))) {
					T::leftChild(T::parent(p), q);
				} else {
					T::rightChild(T::parent(p), q);
				}
			} //else set q as root (p was old root)
			if(q) {
				T::parent(q, T::parent(p));
			}
		}
		
		static void setLeft(Node* const p, Node* const q) {
			if(p) {
				T::leftChild(p, q);
			}
			if(q) {
				T::parent(q, p);
			}
		}
		
		static void setRight(Node* const p, Node* const q) {
			if(p) {
				T::rightChild(p, q);
			}
			if(q) {
				T::parent(q, p);
			}
		}
		
		static Height calcHeight(const Node* const p) {
			return 
				std::max(
					T::leftChild(p)  ? T::height(T::leftChild(p))  : 0, 
					T::rightChild(p) ? T::height(T::rightChild(p)) : 0
				) + 1;
		}
		
		static Size calcVolume(const Node* const p) {
			return 
				(T::leftChild(p)  ? T::volume(T::leftChild(p))  : 0) + 
				(T::rightChild(p) ? T::volume(T::rightChild(p)) : 0) + 1;
		}

		static void augment(Node* const p) {
			T::height(p, calcHeight(p));
			T::volume(p, calcVolume(p));
			T::augment(p);
		}
		
		static void augmentUp(Node* const p) {
			if(p) {
				augment(p);
				augmentUp(T::parent(p));//expects TCO
			}
		}
		
		static void rotateRight(Node* const p) {
			Node* const q = T::leftChild(p);
			T::leftChild(p, T::rightChild(q));
			if(T::leftChild(p)) {
				T::parent(T::leftChild(p), p);
			}
			T::parent(q, T::parent(p));
			if(Node* qp = T::parent(q)) {
				if(T::rightChild(qp) == p) {
					T::rightChild(qp, q);
				} else {
					T::leftChild(qp, q);
				}
			}// else set q as the new root (p was the old root)
			T::rightChild(q, p);
			T::parent(p, q);
			augment(p);
			augment(q);
		}
		
		static void rotateLeft(Node* const p) {
			Node* const q = T::rightChild(p);
			T::rightChild(p, T::leftChild(q));
			if(T::rightChild(p)) {
				T::parent(T::rightChild(p), p);
			}
			T::parent(q, T::parent(p));
			if(Node* qp = T::parent(q)) {
				if(T::leftChild(qp) == p) {
					T::leftChild(qp, q);
				} else {
					T::rightChild(qp, q);
				}
			}// else set q as the new root (p was the old root)
			T::leftChild(q, p);
			T::parent(p, q);
			augment(p);
			augment(q);
		}
		
		static int slope(const Node* const p) {
			return 
				(T::leftChild(p)  ? T::height(T::leftChild(p))  : 0) - 
				(T::rightChild(p) ? T::height(T::rightChild(p)) : 0);
		}
		
		static void balance(Node* const p) {
			if(! p) {
				return;
			}
			augment(p);
			const int bf = slope(p);
			if(bf > 1) {
				if(slope(T::leftChild(p)) == -1) {
					rotateLeft(T::leftChild(p));
				}
				rotateRight(p);
			} else if(bf < -1) {
				if(slope(T::rightChild(p)) == 1) {
					rotateRight(T::rightChild(p));
				}
				rotateLeft(p);
			}
			balance(T::parent(p));//expects TCO
		}
		
		template<class CNode>
		static CNode* findRoot(CNode* p) {
			while(T::parent(p)) {
				p = T::parent(p);
			}
			return p;
		}
		
		template<class CNode>
		static CNode* findHead(CNode* p) {
			while(T::leftChild(p)) {
				p = T::leftChild(p);
			}
			return p;
		}
		
		template<class CNode>
		static CNode* findTail(CNode* p) {
			while(T::rightChild(p)) {
				p = T::rightChild(p);
			}
			return p;
		}
		
		static void insertNodeBefore(Node* const p, Node* const q) {
			T::leftChild(q, nullptr);
			T::rightChild(q, nullptr);
			T::height(q, 1);
			T::volume(q, 1);
			if(T::leftChild(p)) {
				Node* const r = findTail(T::leftChild(p));
				setRight(r, q);
				balance(r);
			} else {
				setLeft(p, q);
				balance(p);
			}
		}
		
		static void insertNodeAfter(Node* const p, Node* const q) {
			T::leftChild(q, nullptr);
			T::rightChild(q, nullptr);
			T::height(q, 1);
			T::volume(q, 1);
			if(T::rightChild(p)) {
				Node* const r = findHead(T::rightChild(p));
				setLeft(r, q);
				balance(r);
			} else {
				setRight(p, q);
				balance(p);
			}
		}
		
		static void removeNode(Node* const p) {
			if(! T::leftChild(p)) {
				Node* const b = T::parent(p);
				transplant(p, T::rightChild(p));
				balance(b);
			} else if(! T::rightChild(p)) {
				Node* const b = T::parent(p);
				transplant(p, T::leftChild(p));
				balance(b);
			} else {
				Node* const q = findHead(T::rightChild(p));
				Node* b = q;
				if(T::parent(q) != p) {
					b = T::parent(q);
					setLeft(T::parent(q), T::rightChild(q));
					setRight(q, T::rightChild(p));
				}
				transplant(p, q);
				setLeft(q, T::leftChild(p));
				balance(b);
			}
			T::parent(p, nullptr);
			T::leftChild(p, nullptr);
			T::rightChild(p, nullptr);
			T::height(p, 1);
			T::volume(p, 1);
		}
		
		static void embed(Node* const p, Node* const q) {
			const Height h = T::height(p);
			Node* const n = findTail(p);
			Node* b = n;
			if(T::parent(n)) {
				b = T::parent(n);
				transplant(n, T::leftChild(n));
				augmentUp(T::parent(n));
			}
			Node* m = q;
			while(h < T::height(m) && T::leftChild(m)) {
				m = T::leftChild(m);
			}
			setLeft(T::parent(m), n);
			setRight(n, m);
			if(n != p) {
				setLeft(n, p);
			}
			balance(b);
		}
		
		static void embrace(Node* const p, Node* const q) {
			const Height h = T::height(q);
			Node* const n = findHead(q);
			Node* b = n;
			if(T::parent(n)) {
				b = T::parent(n);
				transplant(n, T::rightChild(n));
				augmentUp(b);
			}
			Node* m = p;
			while(h < T::height(m) && T::rightChild(m)) {
				m = T::rightChild(m);
			}
			setRight(T::parent(m), n);
			setLeft(n, m);
			if(n != q) {
				setRight(n, q);
			}
			balance(b);
		}
		
		static void join(Node* const p, Node* const q) {
			Node* const rp = findRoot(p);
			Node* const rq = findRoot(q);
			if(T::height(rp) < T::height(rq)) {
				embed(rp, rq);
			} else {
				embrace(rp, rq);
			}
		}
		
		static void makeRoot(Node* const p) {
			if(T::parent(p)) {
				if(T::leftChild(T::parent(p)) == p) {
					rotateRight(T::parent(p));
				} else {
					rotateLeft(T::parent(p));
				}
				makeRoot(p);//expects TCO
			}
		}
		
		static void augmentHeightUp(Node* const p) {
			if(p) {
				Height h = calcHeight(p);
				if(T::height(p) != h) {
					T::height(p, h);
					augmentHeightUp(T::parent(p));//expects TCO
				}
			}
		}
		
		static void balanceDown(Node* const p) {
			if(! p) {
				return;
			}
			augment(p);
			const int bf = slope(p);
			if(bf > 1) {
				if(slope(T::leftChild(p)) == -1) {
					rotateLeft(T::leftChild(p));
				}
				rotateRight(p);
				balanceDown(p);//expects TCO
			} else if(bf < -1) {
				if(slope(T::rightChild(p)) == 1) {
					rotateRight(T::rightChild(p));
				}
				rotateLeft(p);
				balanceDown(p);//expects TCO
			} else {
				balanceDown(T::parent(p));//expects TCO
			}
		}
		
		static void splitBefore(Node* const p) {
			makeRoot(p);
			Node* const q = T::leftChild(p);
			if(q) {
				T::leftChild(p, nullptr);
				T::parent(q, nullptr);
				balanceDown(findTail(q));
			}
			balanceDown(p);
		}
		
		static void splitAfter(Node* const p) {
			makeRoot(p);
			Node* const q = T::rightChild(p);
			if(q) {
				T::rightChild(p, nullptr);
				T::parent(q, nullptr);
				balanceDown(findHead(q));
			}
			balanceDown(p);
		}
		
		static Size size(const Node* const p) {
			return p ? T::volume(p) : 0;
		}
		
		template<class CNode>
		static Node* aboveRightTop(CNode* p) {
			if(! T::parent(p)) {
				return nullptr;
			} else if(p == T::leftChild(T::parent(p))) {
				return T::parent(p);
			} else {
				return aboveRightTop(T::parent(p));//expects TCO
			}
		}
		
		template<class CNode>
		static optional<CNode> next(CNode* p) {
			if(T::rightChild(p)) {
				p = findHead(T::rightChild(p));
			} else {
				p = aboveRightTop(p);
			}
			return p ? 
				optional<CNode>(*p) : 
				optional<CNode>();
		}
		
		template<class CNode>
		static Node* aboveLeftTop(CNode* p) {
			if(! T::parent(p)) {
				return nullptr;
			} else if(p == T::rightChild(T::parent(p))) {
				return T::parent(p);
			} else {
				return aboveLeftTop(T::parent(p));//expects TCO
			}
		}
		
		template<class CNode>
		static optional<CNode> previous(CNode* p) {
			if(T::leftChild(p)) {
				p = findTail(T::leftChild(p));
			} else {
				p = aboveLeftTop(p);
			}
			return p ? 
				optional<CNode>(*p) : 
				optional<CNode>();
		}
		
	
	#ifdef CNRC_DEBUG_AVL_SEQUENCE // helpers for debugging
		
	public:
			
		struct Debug {
			/**
			 * Checks AVL invariants of the tree rooted at 'p'.
			 * Compile with -DDEBUG_CNRC flag to use this.
			 */
			template<class ShowNode>
			static bool checkSanity(
					const Node& p, 
					ShowNode showNode
			) {
				return checkSanity(&p, T::parent(&p), showNode);
			}
			
			template<class ShowNode>
			static bool checkSanity(
				const Node* const p, 
				const Node* const parent, 
				ShowNode showNode
			) {
				if(! p) {
					return true;
				}
				if(T::parent(p) != parent) {
					std::cerr << 
						"AVLSequence sanity check: " << 
						"invalid parent pointer " <<
						"at node " << showNode(*p) << std::endl;
					return false;
				}
				if(T::height(p) != calcHeight(p)) {
					std::cerr << 
						"AVLSequence sanity check: " << 
						"invalid height value " <<
						"at node " << showNode(*p) << std::endl;
					return false;
				}
				const int b = slope(p);
				if(b < -1 || b > 1) {
					std::cerr << 
						"AVLSequence sanity check: " << 
						"violation of the balance condition " <<
						"at node " << showNode(*p) << std::endl;
					return false;
				}
				if(T::volume(p) != calcVolume(p)) {
					std::cerr << 
						"AVLSequence sanity check: " << 
						"invalid branch size value " <<
						"at node: " << showNode(*p) << std::endl;
					return false;
				}
				return 
					checkSanity(T::leftChild(p), p, showNode) && 
					checkSanity(T::rightChild(p), p, showNode);
			}
			
			/**
			 * Prints the structure of the tree rootead at 'p'.
			 * Compile with -DDEBUG_CNRC flag to use this.
			 */
			template<class ShowNode>
			static void printTree(
				const Node& p, 
				ShowNode showNode
			) {
				printTree(&p, "[T]", 0, showNode);
			}
			
			template<class ShowNode>
			static void printTree(
				const Node* const p, 
				const std::string head, 
				const unsigned int indent, 
				ShowNode showNode
			) {
				if(! p) {
					return;
				}
				for(unsigned int i = 0; i < indent; ++i) {
					std::cerr << "\t";
				}
				std::cerr << 
					head << showNode(*p) << ", " <<
					"height: " << T::height(p) << ", " <<
					"branch size: " << T::volume(p) << ", " <<
					"addr: " << p << ", " << 
					"parent addr: " << T::parent(p) << std::endl;
				printTree(T::leftChild(p), "[L]", indent + 1, showNode);
				printTree(T::rightChild(p), "[R]", indent + 1, showNode);
			}
		};
		
	#endif //end CNRC_DEBUG_AVL_SEQUENCE
		
	};
};
	

/**
 * Public interface exported from AVLSequenceModule.
 */
 
template<class NodeTrait> 
using AVLSequenceAlgorithm = AVLSequenceModule::AVLSequenceAlgorithm<NodeTrait>;

template<class T>
using AVLSequenceNodeTrait = AVLSequenceModule::AVLSequenceNodeTrait<T>;

template<class T, class SizeType = unsigned int>
using AVLSequenceNodeMixin = AVLSequenceModule::AVLSequenceNodeMixin<T, SizeType>;

} //end namespace AVLSequence
} //end namespace Cnrc
} //end namespace Snu
