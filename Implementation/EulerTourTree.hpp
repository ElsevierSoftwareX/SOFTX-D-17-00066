/**
 * An implementation of the Euler tour tree described in [J. ACM, 1999, Henzinger and King].
 * This implementation doesn't augment the number of active occurrences.
 */

#pragma once

#include <iterator>
#include <vector>
#include <type_traits>
#include "AVLSequence.hpp"

namespace Snu {
namespace Cnrc {
namespace EulerTourTree {
	
class EulerTourTreeModule {
	
public:
	
	using Size = unsigned int;
	
private:
	
	template<class S, class T>
	using SequenceNodeImpl = AVLSequence::AVLSequenceNodeMixin<S, T>;
	
	template<class T>
	using SequenceAlgorithmImpl = AVLSequence::AVLSequenceAlgorithm<T>;
	
public:
	
	template<class, class> class EulerTourTreeNodeMixin;
	template<class, class> class EulerTourTreeEdgeMixin;
	
private:
	
	template<class Node, class Edge>
	class Occurrence : public SequenceNodeImpl<Occurrence<Node, Edge>, Size> {
		
	public:
		
		Occurrence(Node* const n, bool a)
		: node(n), leftEdge(nullptr), rightEdge(nullptr), isActive(a) {
		}
		
		Occurrence(const Occurrence&) = delete;
		
	public:
		
		Occurrence<Node, Edge>& operator=(const Occurrence<Node, Edge>&) = delete;
		
		Occurrence<Node, Edge>& operator=(Occurrence<Node, Edge>&&) = delete;
		
		bool operator==(const Occurrence<Node, Edge>& other) const {
			return this == &other;
		}
		
		bool operator!=(const Occurrence<Node, Edge>& other) const {
			return this != &other;
		}
		
	public:
		
		Node* const node;
		Edge* leftEdge;
		Edge* rightEdge;
		bool isActive;
	};

public:

	/**
	 * This trait defines an interface for the type `Node` and `Edge` to be stored in an
	 * intrusive Euler tour tree.
	 * If `Node` and `Edge` extend EulerTourTreeNodeMixin<Node, Edge> and EulerTourTreeEdgeMixin<Node, Edge>,
	 * then the default implementation will just work.
	 * If not, a specialization of this template for `Node` and `Edge` must be provided.
	 */
	template<class Node, class Edge>
	class EulerTourTreeTrait {

	private:

		using NodeMixin = EulerTourTreeNodeMixin<Node, Edge>;
		using EdgeMixin = EulerTourTreeEdgeMixin<Node, Edge>;

		static_assert(
			std::is_base_of<NodeMixin, Node>::value, 
			"The default implementation of EulerTourTreeTrait<S, T> requires "
			"S extending EulerTourTreeNodeMixin<S, T>. "
			"If you don't want this, you need to provide a specialization of EulerTourTreeTrait for S, T."
		);

		static_assert(
			std::is_base_of<EdgeMixin, Edge>::value, 
			"The default implementation of EulerTourTreeTrait<S, T> requires "
			"T extending EulerTourTreeEdgeMixin<S, T>. "
			"If you don't want this, you need to provide a specialization of EulerTourTreeTrait for S, T."
		);

	public:

		static Occurrence<Node, Edge>& activeOccurrence(Node& n) {
			return *static_cast<NodeMixin&>(n).activeOccurrence;
		}

		static const Occurrence<Node, Edge>& activeOccurrence(const Node& n) {
			return *static_cast<const NodeMixin&>(n).activeOccurrence;
		}

		static void activeOccurrence(Node& n, Occurrence<Node, Edge>& occ) {
			static_cast<NodeMixin&>(n).activeOccurrence = &occ;
		}

		static void invalidate(Node& n) {
			static_cast<NodeMixin&>(n).activeOccurrence = nullptr;
		}

		static Occurrence<Node, Edge>& occurrence1(Edge& e) {
			return *static_cast<EdgeMixin&>(e).occurrence1;
		}

		static void occurrence1(Edge& e, Occurrence<Node, Edge>& occ) {
			static_cast<EdgeMixin&>(e).occurrence1 = &occ;
		}

		static Occurrence<Node, Edge>& occurrence2(Edge& e) {
			return *static_cast<EdgeMixin&>(e).occurrence2;
		}

		static void occurrence2(Edge& e, Occurrence<Node, Edge>& occ) {
			static_cast<EdgeMixin&>(e).occurrence2 = &occ;
		}

		static Occurrence<Node, Edge>& occurrence3(Edge& e) {
			return *static_cast<EdgeMixin&>(e).occurrence3;
		}

		static void occurrence3(Edge& e, Occurrence<Node, Edge>& occ) {
			static_cast<EdgeMixin&>(e).occurrence3 = &occ;
		}

		static Occurrence<Node, Edge>& occurrence4(Edge& e) {
			return *static_cast<EdgeMixin&>(e).occurrence4;
		}

		static void occurrence4(Edge& e, Occurrence<Node, Edge>& occ) {
			static_cast<EdgeMixin&>(e).occurrence4 = &occ;
		}

		static void invalidate(Edge& e) {
			static_cast<EdgeMixin&>(e).occurrence1 = nullptr;
			static_cast<EdgeMixin&>(e).occurrence2 = nullptr;
			static_cast<EdgeMixin&>(e).occurrence3 = nullptr;
			static_cast<EdgeMixin&>(e).occurrence4 = nullptr;
		}
	};

	/**
	 * This class works as a mixin by CRTP for instances of `Node` to be stored in 
	 * an Euler tour tree connected by the instances of `Edge`.
	 */
	template<class Node, class Edge>
	class EulerTourTreeNodeMixin {

	private:
		using Occ = Occurrence<Node, Edge>;
		
	public:
		
		EulerTourTreeNodeMixin()
		: activeOccurrence(new Occ(static_cast<Node*>(this), true)) {
		}
		
		~EulerTourTreeNodeMixin() {
			using SequenceAlgorithm = SequenceAlgorithmImpl<Occ>;
			if(activeOccurrence) {
				Occ& r = SequenceAlgorithm::findRoot(*activeOccurrence);
				std::vector<Occ*> occs(SequenceAlgorithm::size(r));
				for(Occ& o : SequenceAlgorithm::containerView(r)) {
					occs.push_back(&o);
					if(o.isActive) {
						EulerTourTreeTrait<Node, Edge>::invalidate(*o.node);
					}
				}
				for(Occ* o : occs) {
					delete o;
				}
			}
		}
		
		EulerTourTreeNodeMixin(const EulerTourTreeNodeMixin<Node, Edge>&) = delete;
		
		EulerTourTreeNodeMixin(EulerTourTreeNodeMixin<Node, Edge>&&) = delete;
		
		EulerTourTreeNodeMixin<Node, Edge>& operator=(const EulerTourTreeNodeMixin<Node, Edge>&) = delete;
		
		EulerTourTreeNodeMixin<Node, Edge>& operator=(EulerTourTreeNodeMixin<Node, Edge>&&) = delete;
		
	private:
		
		Occ* activeOccurrence;
			
		template<class, class> friend class EulerTourTreeTrait;	
	};

	/**
	 * This class works as a mixin by CRTP for instances of `Edge` to work as edges 
	 * in an Euler tour tree connecting the instances of `Node`.
	 */
	template<class Node, class Edge>
	class EulerTourTreeEdgeMixin {

	private:
		using Occ = Occurrence<Node, Edge>;
		
	public:
		
		EulerTourTreeEdgeMixin() {
		}
		
		EulerTourTreeEdgeMixin(const EulerTourTreeEdgeMixin<Node, Edge>&) = delete;
		
		EulerTourTreeEdgeMixin(EulerTourTreeEdgeMixin<Node, Edge>&& other) = delete;
		
		EulerTourTreeEdgeMixin<Node, Edge>& operator=(const EulerTourTreeEdgeMixin<Node, Edge>&) = delete;
		
		EulerTourTreeEdgeMixin<Node, Edge>& operator=(EulerTourTreeEdgeMixin<Node, Edge>&& other) = delete;
		
	private:
		
		Occ* occurrence1;
		Occ* occurrence2;
		Occ* occurrence3;
		Occ* occurrence4;
			
		template<class, class> friend class EulerTourTreeTrait;
	};

	/**
	 * You can just use EulerTourTreeNode and EulerTourTreeEdge
	 * instead of writing classes extending EulerTourTreeNodeMixin and EulerTourTreeEdgeMixin,
	 * if you don't need to equip additional data or functionality.
	 */
	struct EulerTourTreeNode;
	struct EulerTourTreeEdge;

	struct EulerTourTreeNode : public EulerTourTreeNodeMixin<EulerTourTreeNode, EulerTourTreeEdge> {
	};

	struct EulerTourTreeEdge : public EulerTourTreeEdgeMixin<EulerTourTreeNode, EulerTourTreeEdge> {
	};
	

	/*
	 * This provides operations on a Euler tour tree of `Node` and `Edge`.
	 * `Node` and `Edge` must implement the interface defined by the template class EulerTourTreeTrait.
	 * Or you can provide an equivalent trait implementation by specifying the third template paramater.
	 */
	template<class Node, class Edge, class ETTTrait = EulerTourTreeTrait<Node, Edge>>
	class EulerTourTreeAlgorithm {

	public:
		
		using Size = EulerTourTreeModule::Size;
	
	private:
		
		using Occ = Occurrence<Node, Edge>;
		using SequenceAlgorithm = SequenceAlgorithmImpl<Occ>;
		
		
	private:
	
		template<class CEdge, class NodeSequenceItr>
		class EdgeIteratorImpl : public std::iterator<std::forward_iterator_tag, CEdge> {
			
		private:
			
			using This = EdgeIteratorImpl<CEdge, NodeSequenceItr>;
			
		public:
			
			EdgeIteratorImpl(NodeSequenceItr i, NodeSequenceItr e)
			: backend(i), end(e) {
				if(i != e) {
					operator++();
				}
			}
			
		public:
			
			This& operator++() {
				++backend;
				while(backend != end &&
					  (ETTTrait::occurrence2(*backend->leftEdge) != *backend)
				) {
					++backend;
				}
				return *this;
			}
			
			This operator++(int) {
				This tmp(*this);
				operator++();
				return tmp;
			}
			
			bool operator==(const This& other) const {
				return backend == other.backend;
			}
			
			bool operator!=(const This& other) const {
				return backend != other.backend;
			}
			
			CEdge& operator*() const {
				return static_cast<Edge&>(*(backend->leftEdge));
			}
			
			CEdge* operator->() const {
				return static_cast<Edge*>(backend->leftEdge);
			}
			
		private:
			NodeSequenceItr backend;
			NodeSequenceItr end;
		};
		
		
		template<class CNode, class COcc>
		class NodeIteratorImpl : public std::iterator<std::forward_iterator_tag, CNode> {
		
		private:
			
			using This = NodeIteratorImpl<CNode, COcc>;
			
		public:
			
			NodeIteratorImpl() 
			: p(nullptr), q(nullptr) {
			}
			
			NodeIteratorImpl(COcc* op, COcc* oq)
			: p(op), q(oq) {
				if(! q) {
					while(p && ! p->isActive) {
						stepForward();
					}
				}
			}
			
		public:
			
			This& operator++() {
				stepForward();
				while(p && ! p->isActive) {
					stepForward();
				}
				return *this;
			}
			
			This operator++(int) {
				This tmp(*this);
				operator++();
				return tmp;
			}
			
			bool operator==(const This& other) const {
				return p == other.p && q == other.q;
			}
			
			bool operator!=(const This& other) const {
				return p != other.p || q != other.q;
			}
			
			CNode& operator*() const {
				return *p->node;
			}
			
			CNode* operator->() const {
				return p->node;
			}
			
		private:
			void stepForward() {
				if(! p) {
					p = q;
					q = nullptr;
				} else if(auto r = SequenceAlgorithm::next(*p)) {
					p = &(r->get());
				} else {
					q = p;
					p = nullptr;
				}
			}
		
		private:
			
			COcc* p;
			COcc* q;
		};
		
	public:
		
		class NodeContainerView {
			
		public:
			
			using iterator = NodeIteratorImpl<Node, Occ>;
			using const_iterator = NodeIteratorImpl<const Node, const Occ>;
			
		public:
			
			explicit NodeContainerView(Occ& root)
			: root(root) {
			}
			
		public:
			
			iterator begin() {
				return iterator(&SequenceAlgorithm::findHead(root), nullptr);
			}
			
			iterator end() {
				return iterator(nullptr, &SequenceAlgorithm::findTail(root));
			}
			
			const_iterator begin() const {
				return const_iterator(&SequenceAlgorithm::findHead(root), nullptr);
			}
			
			const_iterator end() const {
				return const_iterator(nullptr, &SequenceAlgorithm::findTail(root));
			}
			
			const_iterator cbegin() const {
				return const_iterator(&SequenceAlgorithm::findHead(root), nullptr);
			}
			
			const_iterator cend() const {
				return const_iterator(nullptr, &SequenceAlgorithm::findTail(root));
			}
			
		private:
			
			Occ& root;
		};
		
		class ConstNodeContainerView {
			
		public:
			
			using const_iterator = NodeIteratorImpl<const Node, const Occ>;
			
		public:
			
			explicit ConstNodeContainerView(const Occ& root)
			: root(root) {
			}
			
		public:
			
			const_iterator begin() const {
				return const_iterator(&SequenceAlgorithm::findHead(root), nullptr);
			}
			
			const_iterator end() const {
				return const_iterator(nullptr, &SequenceAlgorithm::findTail(root));
			}
			
			const_iterator cbegin() const {
				return const_iterator(&SequenceAlgorithm::findHead(root), nullptr);
			}
			
			const_iterator cend() const {
				return const_iterator(nullptr, &SequenceAlgorithm::findTail(root));
			}
			
		private:
			
			const Occ& root;
		};
		
		class EdgeContainerView {
			
		private:
			
			using SequenceContainerView = typename SequenceAlgorithm::ContainerView;
			using SequenceNodeIterator = typename SequenceContainerView::iterator;
			using ConstSequenceNodeIterator = typename SequenceContainerView::const_iterator;
			
		public:
			
			using iterator = EdgeIteratorImpl<Edge, SequenceNodeIterator>;
			using const_iterator = EdgeIteratorImpl<const Edge, ConstSequenceNodeIterator>;
			
		public:
			
			explicit EdgeContainerView(SequenceContainerView b)
			: backend(b) {
			}
			
		public:
			
			iterator begin() {
				return iterator(backend.begin(), backend.end());
			}
			
			iterator end() {
				return iterator(backend.end(), backend.end());
			}
			
			const_iterator begin() const {
				return const_iterator(backend.cbegin(), backend.cend());
			}
			
			const_iterator end() const {
				return const_iterator(backend.cend(), backend.cend());
			}
			
			const_iterator cbegin() const {
				return const_iterator(backend.cbegin(), backend.cend());
			}
			
			const_iterator cend() const {
				return const_iterator(backend.cend(), backend.cend());
			}
			
			
		private:
			
			SequenceContainerView backend;
		};
		
		class ConstEdgeContainerView {
			
		private:
			
			using ConstSequenceContainerView = typename SequenceAlgorithm::ConstContainerView;
			using ConstSequenceNodeIterator = typename ConstSequenceContainerView::const_iterator;
			
		public:
			
			using const_iterator = EdgeIteratorImpl<const Edge, ConstSequenceNodeIterator>;
			
		public:
			
			explicit ConstEdgeContainerView(ConstSequenceContainerView b)
			: backend(b) {
			}
			
		public:
			
			const_iterator begin() const {
				return const_iterator(backend.begin(), backend.end());
			}
			
			const_iterator end() const {
				return const_iterator(backend.end(), backend.end());
			}
			
			const_iterator cbegin() const {
				return const_iterator(backend.cbegin(), backend.end());
			}
			
			const_iterator cend() const {
				return const_iterator(backend.cend(), backend.end());
			}
			
		private:
			
			ConstSequenceContainerView backend;
		};
		
	public:
		
		/**
		 * Returns a NodeContainerView of the tree of `n`.
		 * It has an interface for STL compatible iterators.
		 * The iterators visit each node once and only once.
		 */
		static NodeContainerView nodeContainerView(Node& n) {
			return NodeContainerView(
				SequenceAlgorithm::findRoot(
					ETTTrait::activeOccurrence(n)
				)
			);
		}
		
		/**
		 * Returns a ConstNodeContainerView of the tree of `n`.
		 * It has an interface for STL compatible iterators.
		 * The iterators visit each node once and only once.
		 */
		static ConstNodeContainerView nodeContainerView(const Node& n) {
			return ConstNodeContainerView(
				SequenceAlgorithm::findRoot(
					ETTTrait::activeOccurrence(n)
				)
			);
		}
		
		/**
		 * Returns a EdgeContainerView of the tree of `n`.
		 * It has an interface for STL compatible iterators.
		 * The iterators visit each node once and only once.
		 */
		static EdgeContainerView edgeContainerView(Node& n) {
			return EdgeContainerView(
				SequenceAlgorithm::containerView(
					SequenceAlgorithm::findRoot(
						ETTTrait::activeOccurrence(n)
					)
				)
			);
		}
		
		/**
		 * Returns a ConstEdgeContainerView of the tree of `n`.
		 * It has an interface for STL compatible iterators.
		 * The iterators visit each node once and only once.
		 */
		static ConstEdgeContainerView edgeContainerView(const Node& n) {
			return ConstEdgeContainerView(
				SequenceAlgorithm::containerView(
					SequenceAlgorithm::findRoot(
						ETTTrait::activeOccurrence(n)
					)
				)
			);
		}
		
		/**
		 * Connects two nodes `n1` and `n2` by the edge `e`.
		 */
		static void connect(Node& n1, Node& n2, Edge& e) {
			Occ& o1h = makeHead(n1);
			Occ& o1r = SequenceAlgorithm::findRoot(o1h);
			Occ& o1t = SequenceAlgorithm::findTail(o1r);
			Occ& o2h = makeHead(n2);
			Occ& o2r = SequenceAlgorithm::findRoot(o2h);
			Occ& o2t = SequenceAlgorithm::findTail(o2r);
			SequenceAlgorithm::join(o1t, o2h);
			Occ& ont = *(new Occ(&n1, false));
			SequenceAlgorithm::insertNodeAfter(o2t, ont);
			o1t.rightEdge = &e;
			o2h.leftEdge  = &e;
			ETTTrait::occurrence1(e, o1t);
			ETTTrait::occurrence2(e, o2h);
			ont.leftEdge  = &e;
			o2t.rightEdge = &e;
			ETTTrait::occurrence3(e, o2t);
			ETTTrait::occurrence4(e, ont);
		}
		
		/**
		 * Disconnects the edge `e`.
		 */
		static void disconnect(Edge& e) {
			SequenceAlgorithm::splitAfter(ETTTrait::occurrence1(e));
			SequenceAlgorithm::splitAfter(ETTTrait::occurrence3(e));
			if( SequenceAlgorithm::findRoot(ETTTrait::occurrence1(e)) == 
				SequenceAlgorithm::findRoot(ETTTrait::occurrence4(e))
			) {
				join(ETTTrait::occurrence3(e), ETTTrait::occurrence2(e));
				ETTTrait::occurrence4(e).leftEdge  = nullptr;
				ETTTrait::occurrence1(e).rightEdge = nullptr;
			} else {
				join(ETTTrait::occurrence1(e), ETTTrait::occurrence4(e));
				ETTTrait::occurrence2(e).leftEdge  = nullptr;
				ETTTrait::occurrence3(e).rightEdge = nullptr;
			}
			ETTTrait::invalidate(e);
		}
		
		/**
		 * Returns only when `n1` and `n2` are in the same tree.
		 */
		static bool hasPath(const Node& n1, const Node& n2) {
			return 
				SequenceAlgorithm::findRoot(ETTTrait::activeOccurrence(n1)) ==
				SequenceAlgorithm::findRoot(ETTTrait::activeOccurrence(n2));
		}
		
		/**
		 * Returns whether of not `n` is the representative of the tree.
		 * Each tree has a representative.
		 * The representative can be changed by any operations modifying the tree.
		 */
		static bool isClusterRep(const Node& n) {
			return &n == &findClusterRep(n);
		}
		
		/**
		 * Returns the representative of the tree of `n`.
		 */
		static Node& findClusterRep(Node& n) {
			return *SequenceAlgorithm::findRoot(ETTTrait::activeOccurrence(n)).node;
		}
		
		/**
		 * Returns the representative of the tree of `n`.
		 */
		static const Node& findClusterRep(const Node& n) {
			return *SequenceAlgorithm::findRoot(ETTTrait::activeOccurrence(n)).node;
		}
		
		/**
		 * Returns the number of distinctive nodes in the tree.
		 */
		static Size clusterSize(const Node& n) {
			const Occ& r = SequenceAlgorithm::findRoot(ETTTrait::activeOccurrence(n));
			return (SequenceAlgorithm::size(r) + 1) / 2;
		}
		
	private:
		
		static void join(Occ& p, Occ& q) {
			if(auto pp = SequenceAlgorithm::previous(p)) {
				SequenceAlgorithm::removeNode(p);
				SequenceAlgorithm::join(*pp, q);
				putOccurrencOnEdge(*pp, q);
			}
			if(p.isActive) {
				ETTTrait::activeOccurrence(*p.node, q);
				q.isActive = true;
			}
			delete &p;
		};
		
		static Occ& makeHead(Node& n) {
			Occ& newHead = ETTTrait::activeOccurrence(n);
			Occ& oldRoot = SequenceAlgorithm::findRoot(newHead);
			Occ& oldHead = SequenceAlgorithm::findHead(oldRoot);
			if(oldHead.node == &n) {
				return oldHead;
			}
			Occ& oldTail = SequenceAlgorithm::findTail(oldRoot);
			Occ& leftOfOldTail = *SequenceAlgorithm::previous(oldTail);
			Occ& leftOfNewTail = *SequenceAlgorithm::previous(newHead);
			Occ& newTail = *(new Occ(&n, false));
			SequenceAlgorithm::splitBefore(newHead);
			SequenceAlgorithm::removeNode(oldTail);
			SequenceAlgorithm::join(leftOfOldTail, oldHead);
			SequenceAlgorithm::insertNodeAfter(leftOfNewTail, newTail);
			newHead.leftEdge = nullptr;
			putOccurrencOnEdge(leftOfOldTail, oldHead);
			putOccurrencOnEdge(leftOfNewTail, newTail);
			if(oldTail.isActive) {
				ETTTrait::activeOccurrence(*oldTail.node, oldHead);
				oldHead.isActive = true;
			}
			delete &oldTail;
			return newHead;
		}
		
		static void putOccurrencOnEdge(Occ& left, Occ& right) {
			Edge& e = *left.rightEdge;
			right.leftEdge = &e;
			if(&ETTTrait::occurrence1(e) == &left) {
				ETTTrait::occurrence2(e, right);
			} else {
				ETTTrait::occurrence4(e, right);
			}
		}
		
	};

};


/**
 * Public interface from EulerTourTreeModule.
 */
template<class Node, class Edge> 
using EulerTourTreeTrait = EulerTourTreeModule::EulerTourTreeTrait<Node, Edge>;

template<class Node, class Edge> 
using EulerTourTreeAlgorithm = EulerTourTreeModule::EulerTourTreeAlgorithm<Node, Edge>;

template<class Node, class Edge>
using EulerTourTreeNodeMixin = EulerTourTreeModule::EulerTourTreeNodeMixin<Node, Edge>;

template<class Node, class Edge>
using EulerTourTreeEdgeMixin = EulerTourTreeModule::EulerTourTreeEdgeMixin<Node, Edge>;

using EulerTourTreeNode = EulerTourTreeModule::EulerTourTreeNode;

using EulerTourTreeEdge = EulerTourTreeModule::EulerTourTreeEdge;

} //end namespace EulerTourTree
} //end namespace Cnrc
} //end namespace Snu
