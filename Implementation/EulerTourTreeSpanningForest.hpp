/**
 * An implementation of the dynamic spanning forest based on Euler tour trees.
 */

#pragma once

#include <set>
#include <iterator>
#include <type_traits>
#include "EulerTourTree.hpp"

namespace Snu {
namespace Cnrc {
namespace EulerTourTreeSpanningForest {

class EulerTourTreeSpanningForestModule {
	
public:
	
	using Size = unsigned int;
	
public:
	
	template<class> class Cluster;
	template<class> class ConstCluster;
	template<class, class> class EulerTourTreeSpanningForestNodeMixin;
	template<class, class> class EulerTourTreeSpanningForestEdgeMixin;
	template<class, class, class> class EulerTourTreeSpanningForestAlgorithm;
	
public:

	/**
	 * This trait defines an interface for the type `Node` and `Edge` to be stored in an
	 * spanning forest based on the Euler tour tree.
	 * If `Node` and `Edge` extend EulerTourTreeSpanningForestNodeMixin<Node, Edge> and EulerTourTreeSpanningForestEdgeMixin<Node, Edge>
	 * respectively, then the default implementation will just work.
	 * If not, a specialization of this template for `Node` and `Edge` must be provided.
	 */
    template<class Node, class Edge>
    class EulerTourTreeSpanningForestTrait {

	private:
		using ETTFNodeMixin = EulerTourTreeSpanningForestNodeMixin<Node, Edge>;
		using ETTFEdgeMixin = EulerTourTreeSpanningForestEdgeMixin<Node, Edge>;

		static_assert(
			std::is_base_of<ETTFNodeMixin, Node>::value, 
			"The default implementation of EulerTourTreeSpanningForestTrait<S, T> requires "
			"S extending EulerTourTreeSpanningForestNodeMixin<S, T>. "
			"If you don't want this, you need to provide a specialization of EulerTourTreeSpanningForestTrait for S, T."
		);

		static_assert(
			std::is_base_of<ETTFEdgeMixin, Edge>::value, 
			"The default implementation of EulerTourTreeSpanningForestTrait<S, T> requires "
			"T extending EulerTourTreeSpanningForestEdgeMixin<S, T>. "
			"If you don't want this, you need to provide a specialization of EulerTourTreeSpanningForestTrait for S, T."
		);

	public:
        
        using EdgeSet = typename ETTFNodeMixin::EdgeSet;
        
        static EdgeSet& edges(Node& n) {
            return static_cast<ETTFNodeMixin&>(n).edges;
        }

        static const EdgeSet& edges(const Node& n) {
            return static_cast<const ETTFNodeMixin&>(n).edges;
        }

        static Node& node1(Edge& e) {
            return *static_cast<ETTFEdgeMixin&>(e).node1;
        }

        static const Node& node1(const Edge& e) {
            return *static_cast<const ETTFEdgeMixin&>(e).node1;
        }

		static void node1(Edge& e, Node& n) {
            static_cast<ETTFEdgeMixin&>(e).node1 = &n;
        }

        static Node& node2(Edge& e) {
            return *static_cast<ETTFEdgeMixin&>(e).node2;
        }

        static const Node& node2(const Edge& e) {
            return *static_cast<const ETTFEdgeMixin&>(e).node2;
        }

		static void node2(Edge& e, Node& n) {
            static_cast<ETTFEdgeMixin&>(e).node2 = &n;
        }

        static bool isTreeEdge(const Edge& e) {
            return static_cast<const ETTFEdgeMixin&>(e).isTreeEdge;
        }

		static void setTreeEdge(Edge& e, bool b) {
            static_cast<ETTFEdgeMixin&>(e).isTreeEdge = b;
        }

		static bool isValid(const Edge& e) {
			return static_cast<const ETTFEdgeMixin&>(e).node1 != nullptr &&
			       static_cast<const ETTFEdgeMixin&>(e).node2 != nullptr;
		}

		static void invalidate(Edge& e) {
			static_cast<ETTFEdgeMixin&>(e).node1 = nullptr;
			static_cast<ETTFEdgeMixin&>(e).node2 = nullptr;
		}
    };

public:

	/**
	 * This class works as a mixin by CRTP for instances of `Node` to be stored in 
	 * an spanning forest connected by the instances of `Edge`.
	 */
	template<class Node, class Edge>
	class EulerTourTreeSpanningForestNodeMixin : public EulerTourTree::EulerTourTreeNodeMixin<Node, Edge> {
		
	public:
		using EdgeSet = std::set<Edge*>;
	
	public:
		
		EulerTourTreeSpanningForestNodeMixin()
		: edges() {
		}
		
		EulerTourTreeSpanningForestNodeMixin(const EulerTourTreeSpanningForestNodeMixin<Node, Edge>&) = delete;
		
		EulerTourTreeSpanningForestNodeMixin(EulerTourTreeSpanningForestNodeMixin<Node, Edge>&&) = delete;
		
	public:
		
		EulerTourTreeSpanningForestNodeMixin& operator=(const EulerTourTreeSpanningForestNodeMixin<Node, Edge>&) = delete;
		
		EulerTourTreeSpanningForestNodeMixin& operator=(EulerTourTreeSpanningForestNodeMixin<Node, Edge>&&) = delete;
		
	private:
		
		EdgeSet edges;
		
        template<class, class> friend class EulerTourTreeSpanningForestTrait;
	};
	
	/**
	 * This class works as a mixin by CRTP for instances of `Edge` to work as edges 
	 * in an spanning forest connecting the instances of `Node`.
	 */
    template<class Node, class Edge>
	class EulerTourTreeSpanningForestEdgeMixin : public EulerTourTree::EulerTourTreeEdgeMixin<Node, Edge> {
		
	public:
		
		EulerTourTreeSpanningForestEdgeMixin()
		: node1(nullptr), node2(nullptr), isTreeEdge(false) {
		}
		
		EulerTourTreeSpanningForestEdgeMixin(const EulerTourTreeSpanningForestEdgeMixin<Node, Edge>&) = delete;
		
		EulerTourTreeSpanningForestEdgeMixin(EulerTourTreeSpanningForestEdgeMixin<Node, Edge>&&) = delete;
		
	private:
		
		EulerTourTreeSpanningForestEdgeMixin& operator=(const EulerTourTreeSpanningForestEdgeMixin<Node, Edge>&) = delete;
		
		EulerTourTreeSpanningForestEdgeMixin& operator=(EulerTourTreeSpanningForestEdgeMixin<Node, Edge>&&) = delete;
		
	private:
		
		Node* node1;
		Node* node2;
		bool isTreeEdge;
		
		template<class, class> friend class EulerTourTreeSpanningForestTrait;
	};

	/**
	 * You can just use EulerTourTreeSpanningForestNode and EulerTourTreeSpanningForestEdge
	 * instead of writing classes extending EulerTourTreeSpanningForestNode and EulerTourTreeSpanningForestEdge,
	 * if you don't need to equip additional data or functionality.
	 */
	class EulerTourTreeSpanningForestNode;
	class EulerTourTreeSpanningForestEdge;

	class EulerTourTreeSpanningForestNode
	: public EulerTourTreeSpanningForestNodeMixin<
		EulerTourTreeSpanningForestNode, 
		EulerTourTreeSpanningForestEdge
	> {
	};

	class EulerTourTreeSpanningForestEdge
	: public EulerTourTreeSpanningForestEdgeMixin<
		EulerTourTreeSpanningForestNode, 
		EulerTourTreeSpanningForestEdge
	> {
	};
	
public:
	
	/*
	 * This provides operations on a spanning forest of `Node` and `Edge`.
	 * `Node` and `Edge` must implement the interface defined by the template class EulerTourTreeSpanningForestTrait.
	 * Or you can provide an equivalent trait implementation by specifying the third template paramater.
	 */
	template<class Node, class Edge, class SFTrait = EulerTourTreeSpanningForestTrait<Node, Edge>>
	class EulerTourTreeSpanningForestAlgorithm {
		
	private:
		
		using ETTAlgorithm = EulerTourTree::EulerTourTreeAlgorithm<Node, Edge>;
		
	public:
		
		class Cluster : public ETTAlgorithm::NodeContainerView {
			
		private:
			
			explicit Cluster(Node& r)
			: ETTAlgorithm::NodeContainerView(ETTAlgorithm::nodeContainerView(r)), rep(r) {
			}
			
		public:
			
			Node& representative() {
				return rep;
			}
			
			const Node& representative() const {
				return rep;
			}
			
			Size size() const {
				return ETTAlgorithm::clusterSize(rep);
			}

			bool contains(const Node& n) const {
				return ETTAlgorithm::hasPath(rep, n);
			}
		
		private:
			Node& rep;
			
			template<class, class, class> friend class EulerTourTreeSpanningForestAlgorithm;
		};
		
		class ConstCluster : public ETTAlgorithm::ConstNodeContainerView {
			
		private:
			
			explicit ConstCluster(const Node& r)
			: ETTAlgorithm::ConstNodeContainerView(ETTAlgorithm::nodeContainerView(r)), rep(r) {
			}
			
		public:	
			
			const Node& representative() const {
				return rep;
			}
			
			Size size() const {
				return ETTAlgorithm::clusterSize(rep);
			}

			bool contains(const Node& n) const {
				return ETTAlgorithm::hasPath(rep, n);
			}
		
		private:
			const Node& rep;
			
			template<class, class, class> friend class EulerTourTreeSpanningForestAlgorithm;
		};
		
	private:
		
		template<class CEdge, class BackendItr>
		class EdgeIteratorImpl : public std::iterator<std::forward_iterator_tag, CEdge> {
			
		private:
			
			using This = EdgeIteratorImpl<CEdge, BackendItr>;
	
		public:
			
			EdgeIteratorImpl() 
			: itr() {
			}
			
			explicit EdgeIteratorImpl(BackendItr itr)
			: itr(itr) {
			}
			
		public:
			
			This& operator++() {
				++itr;
				return *this;
			}
			
			This operator++(int) {
				This tmp(*this);
				operator++();
				return tmp;
			}
			
			bool operator==(const This& other) const {
				return itr == other.itr;
			}
			
			bool operator!=(const This& other) const {
				return itr != other.itr;
			}
			
			CEdge& operator*() const {
				return **itr;
			}
			
			CEdge* operator->() const {
				return *itr;
			}
		
		private:
			
			BackendItr itr;
		};
		
	public:
		
		class Edges {
			
		private:
			
			explicit Edges(Node& n)
			: node(n) {
			}
			
		public:
			
			using iterator = EdgeIteratorImpl<Edge, typename Node::EdgeSet::iterator>;
			using const_iterator = EdgeIteratorImpl<const Edge, typename Node::EdgeSet::const_iterator>;
			
		public:
			
			iterator begin() {
				return iterator(SFTrait::edges(node).begin());
			}
			
			iterator end() {
				return iterator(SFTrait::edges(node).end());
			}
			
			const_iterator begin() const {
				return const_iterator(SFTrait::edges(node).cbegin());
			}
			
			const_iterator end() const {
				return const_iterator(SFTrait::edges(node).cend());
			}
			
			const_iterator cbegin() const {
				return const_iterator(SFTrait::edges(node).cbegin());
			}
			
			const_iterator cend() const {
				return const_iterator(SFTrait::edges(node).cend());
			}
			
		private:
			Node& node;
			
			template<class, class, class> friend class EulerTourTreeSpanningForestAlgorithm;
		};
		
		class ConstEdges {
			
		private:
			
			explicit ConstEdges(const Node& n)
			: node(n) {
			}
			
		public:
			
			using const_iterator = EdgeIteratorImpl<const Edge, typename Node::EdgeSet::const_iterator>;
			
		public:
			
			const_iterator begin() const {
				return const_iterator(SFTrait::edges(node).cbegin());
			}
			
			const_iterator end() const {
				return const_iterator(SFTrait::edges(node).cend());
			}
			
			const_iterator cbegin() const {
				return const_iterator(SFTrait::edges(node).cbegin());
			}
			
			const_iterator cend() const {
				return const_iterator(SFTrait::edges(node).cend());
			}
			
		private:
			const Node& node;
			
			template<class, class, class> friend class EulerTourTreeSpanningForestAlgorithm;
		};
		
		
	public:
		
		/**
		 * Returns whether or not `n1` and `n2` are in the same connected component.
		 */
		static bool hasPath(const Node& n1, const Node& n2) {
			return ETTAlgorithm::hasPath(n1, n2);
		}
		
		/**
		 * Returns a Cluster which represents the connected compnent of `n`.
		 * It has an interface for STL compatible iterators with which
		 * you can visit each node in the component one by one.
		 */
		static Cluster cluster(Node& n) {
			return Cluster(ETTAlgorithm::findClusterRep(n)); 
		}
		
		/**
		 * Returns a ConstCluster which represents the connected compnent of `n`.
		 * It has an interface for STL compatible iterators with which
		 * you can visit each node in the component one by one.
		 */
		static ConstCluster cluster(const Node& n) {
			return ConstCluster(ETTAlgorithm::findClusterRep(n)); 
		}

		/**
		 * Returns a ConstCluster which represents the connected compnent of `n`.
		 * It has an interface for STL compatible iterators with which
		 * you can visit each node in the component one by one.
		 */
		static ConstCluster constCluster(Node& n) {
			return cluster(static_cast<const Node&>(n)); 
		}
		
		/**
		 * Returns the number of nodes in the connected component of `n`.
		 */
		static Size clusterSize(const Node& n) {
			return ETTAlgorithm::clusterSize(ETTAlgorithm::findClusterRep(n));
		}

		/**
		 * Returns whether of not `n` is the representative of its connected component.
		 * Each connected component has a representative.
		 * The representative can be changed by any operations modifying the connected component.
		 */
		static bool isClusterRep(const Node& n) {
			return &ETTAlgorithm::findClusterRep(n) == &n;
		}
		
		/**
		 * Returns the representative of the connected component of `n`.
		 */
		static const Node& findClusterRep(const Node& n) {
			return ETTAlgorithm::findClusterRep(n);
		}
		
		/**
		 * Returns the representative of the connected component of `n`.
		 */
		static Node& findClusterRep(Node& n) {
			return ETTAlgorithm::findClusterRep(n);
		}
		
		/**
		 * Returns an Edges with which you can visit the edges of `n` one by one.
		 * It has an interface for STL compatible iterators. 
		 */
		static Edges edges(Node& n) {
			return Edges(n);
		}
		
		/**
		 * Returns a ConstEdges with which you can visit the edges of `n` one by one.
		 * It has an interface for STL compatible iterators. 
		 */
		static ConstEdges edges(const Node& n) {
			return ConstEdges(n);
		}

		/**
		 * Returns a ConstEdges with which you can visit the edges of `n` one by one.
		 * It has an interface for STL compatible iterators. 
		 */
		static ConstEdges constEdges(Node& n) {
			return edges(static_cast<const Node&>(n));
		}
		
		/**
		 * Connects two nodes `n1` and `n2` by the edge `e`.
		 */
		static bool connect(Node& n1, Node& n2, Edge& e) {
			SFTrait::node1(e, n1);
			SFTrait::node2(e, n2);
			SFTrait::edges(n1).insert(&e);
			SFTrait::edges(n2).insert(&e);
			if(ETTAlgorithm::hasPath(n1, n2)) {
				SFTrait::setTreeEdge(e, false);
				return false;
			} else {
				replaceWith(e);
				return true;
			}
		}
		
		/**
		 * Disconnects the edge `e` and invalidates it.
		 * You can reuse the `e` by connecting any two nodes.
		 * Any other operations on `e` is not allowed until the reuse.
		 */
		static bool disconnect(Edge& e) {
			bool clusterSplitted = false;
			SFTrait::edges(SFTrait::node1(e)).erase(&e);
			SFTrait::edges(SFTrait::node2(e)).erase(&e);
			if(SFTrait::isTreeEdge(e)) {
				ETTAlgorithm::disconnect(e);
				clusterSplitted = ! checkReplacement(e);
			}
			invalidate(e);
			return clusterSplitted;
		}
		
		/**
		 * Returns the first end of the edge `e`.
		 */
		static const Node& node1(const Edge& e) {
			return SFTrait::node1(e);
		}
		
		/**
		 * Returns the first end of the edge `e`.
		 */
		static Node& node1(Edge& e) {
			return SFTrait::node1(e);
		}
		
		/**
		 * Returns the second end of the edge `e`.
		 */
		static const Node& node2(const Edge& e) {
			return SFTrait::node2(e);
		}
		
		/**
		 * Returns the second end of the edge `e`.
		 */
		static Node& node2(Edge& e) {
			return SFTrait::node2(e);
		}

		/**
		 * Returns the validity of the edge `e`.
		 */
		static bool isValid(const Edge& e) {
			return SFTrait::isValid(e);
		}

		/**
		 * Invalidates the edge `e`.
		 */
		static void invalidate(Edge& e) {
			SFTrait::invalidate(e);
		}

	private:
		
		static void replaceWith(Edge& e) {
			ETTAlgorithm::connect(SFTrait::node1(e), SFTrait::node2(e), e);
			SFTrait::setTreeEdge(e, true);
		}
		
		static bool checkReplacement(Edge& e) {
			Node& r1 = ETTAlgorithm::findClusterRep(SFTrait::node1(e));
			Node& r2 = ETTAlgorithm::findClusterRep(SFTrait::node2(e));
			typename ETTAlgorithm::Size sz1 = ETTAlgorithm::clusterSize(r1);
			typename ETTAlgorithm::Size sz2 = ETTAlgorithm::clusterSize(r2);
			if(sz1 < sz2) {
				return checkReplacement(r1, r2);
			} else {
				return checkReplacement(r2, r1);
			}
		}

		static bool checkReplacement(Node& smaller, Node& larger) {
			for(Node& n : ETTAlgorithm::nodeContainerView(smaller)) {
				for(Edge* e : SFTrait::edges(n)) {
					Node& m = &SFTrait::node1(*e) == &n ? SFTrait::node2(*e) : SFTrait::node1(*e);
					if(&ETTAlgorithm::findClusterRep(m) == &larger) {
						replaceWith(*e);
						return true;
					}
				}
			}
			return false;
		}
	};
	
};


/**
 * Public interface from EulerTourTreeSpanningForestModule.
 */

template<class Node, class Edge>
using EulerTourTreeSpanningForestTrait = EulerTourTreeSpanningForestModule::EulerTourTreeSpanningForestTrait<Node, Edge>;

template<class Node, class Edge>
using EulerTourTreeSpanningForestAlgorithm = EulerTourTreeSpanningForestModule::EulerTourTreeSpanningForestAlgorithm<Node, Edge>;

template<class Node, class Edge>
using EulerTourTreeSpanningForestNodeMixin = EulerTourTreeSpanningForestModule::EulerTourTreeSpanningForestNodeMixin<Node, Edge>;

template<class Node, class Edge>
using EulerTourTreeSpanningForestEdgeMixin = EulerTourTreeSpanningForestModule::EulerTourTreeSpanningForestEdgeMixin<Node, Edge>;

using EulerTourTreeSpanningForestNode = EulerTourTreeSpanningForestModule::EulerTourTreeSpanningForestNode;

using EulerTourTreeSpanningForestEdge = EulerTourTreeSpanningForestModule::EulerTourTreeSpanningForestEdge;


} //end namespace EulerTourTreeSpanningForest
} //end namespace Cnrc
} //end namespace Snu
