/**
 * An implementation of the fully dynamic graph connectivity algorithm 
 * by Holm, de Lichetenberg and Thorup (J. ACM, 2001).
 */

#include <set>
#include <vector>
#include <iterator>
#include <type_traits>
#include <utility>
#include "EulerTourTree.hpp"

namespace Snu {
namespace Cnrc {
namespace HDTSpanningForest {

class HDTSpanningForestModule {
	
public:
	
	using Level = unsigned int;
	using Size = unsigned int;

private:
	template<class, class> class LevelNode;
	template<class, class> class LevelEdge;
	
public:
	
	template<class> class Cluster;
	template<class> class ConstCluster;
	template<class, class> class HDTSpanningForestNodeMixin;
	template<class, class> class HDTSpanningForestEdgeMixin;
	template<class, class, class> class HDTSpanningForestAlgorithm;
	
private:

	template<class Node, class Edge>
	class LevelNode : public EulerTourTree::EulerTourTreeNodeMixin<LevelNode<Node, Edge>, LevelEdge<Node, Edge>> {
		
	public:

		using EdgeSet = std::set<Edge*>;
		
		explicit LevelNode(Node& n)
		: edges(), boxNode(n) {
		}
		
		bool operator==(LevelNode& other) const {
			return this == &other;
		}
		
		bool operator!=(LevelNode& other) const {
			return this != &other;
		}
		
		EdgeSet edges;
		Node& boxNode;
	};
	
	template<class Node, class Edge>
	class LevelEdge : public EulerTourTree::EulerTourTreeEdgeMixin<LevelNode<Node, Edge>, LevelEdge<Node, Edge>> {
		
	public:
	
		explicit LevelEdge(Edge& e)
		: boxEdge(e) {
		}
		
		Edge& boxEdge;
	};

	/**
	 * This trait defines an interface for the type `Node` and `Edge` to be stored in an
	 * spanning forest based on HDT algorithm.
	 * If `Node` and `Edge` extend HDTSpanningForestNodeMixin<Node, Edge> and HDTSpanningForestEdgeMixin<Node, Edge>
	 * respectively, then the default implementation will just work.
	 * If not, a specialization of this template for `Node` and `Edge` must be provided.
	 */
	template<class Node, class Edge>
    class HDTSpanningForestTrait {

	private:
		using HDTNodeMixin = HDTSpanningForestNodeMixin<Node, Edge>;
		using HDTEdgeMixin = HDTSpanningForestEdgeMixin<Node, Edge>;

		static_assert(
			std::is_base_of<HDTNodeMixin, Node>::value, 
			"The default implementation of HDTSpanningForestTrait<S, T> requires "
			"S extending HDTSpanningForestNodeMixin<S, T>. "
			"If you don't want this, you need to provide a specialization of HDTSpanningForestTrait for S, T."
		);

		static_assert(
			std::is_base_of<HDTEdgeMixin, Edge>::value, 
			"The default implementation of HDTSpanningForestTrait<S, T> requires "
			"T extending HDTSpanningForestEdgeMixin<S, T>. "
			"If you don't want this, you need to provide a specialization of HDTSpanningForestTrait for S, T."
		);

	public:
        
        static std::vector<LevelNode<Node, Edge>*>& levelNodes(Node& n) {
			return static_cast<HDTNodeMixin&>(n).levelNodes;
		}

		static const std::vector<LevelNode<Node, Edge>*>& levelNodes(const Node& n) {
			return static_cast<const HDTNodeMixin&>(n).levelNodes;
		}

		static std::vector<LevelEdge<Node, Edge>*>& levelEdges(Edge& e) {
			return static_cast<HDTEdgeMixin&>(e).levelEdges;
		}

		static const std::vector<LevelEdge<Node, Edge>*>& levelEdges(const Edge& e) {
			return static_cast<const HDTEdgeMixin&>(e).levelEdges;
		}

        static Node& node1(Edge& e) {
            return *static_cast<HDTEdgeMixin&>(e).node1;
        }

        static const Node& node1(const Edge& e) {
            return *static_cast<const HDTEdgeMixin&>(e).node1;
        }

		static void node1(Edge& e, Node& n) {
            static_cast<HDTEdgeMixin&>(e).node1 = &n;
        }

        static Node& node2(Edge& e) {
            return *static_cast<HDTEdgeMixin&>(e).node2;
        }

        static const Node& node2(const Edge& e) {
            return *static_cast<const HDTEdgeMixin&>(e).node2;
        }

		static void node2(Edge& e, Node& n) {
            static_cast<HDTEdgeMixin&>(e).node2 = &n;
        }

        static bool isTreeEdge(const Edge& e) {
            return static_cast<const HDTEdgeMixin&>(e).isTreeEdge;
        }

		static void setTreeEdge(Edge& e, bool b) {
            static_cast<HDTEdgeMixin&>(e).isTreeEdge = b;
        }

		static bool isValid(const Edge& e) {
			return static_cast<const HDTEdgeMixin&>(e).node1 != nullptr &&
			       static_cast<const HDTEdgeMixin&>(e).node2 != nullptr;
		}

		static void invalidate(Edge& e) {
			static_cast<HDTEdgeMixin&>(e).node1 = nullptr;
			static_cast<HDTEdgeMixin&>(e).node2 = nullptr;
		}

		static Level level(const Edge& e) {
			return static_cast<const HDTEdgeMixin&>(e).level;
		}

		static void level(Edge& e, Level l) {
			static_cast<HDTEdgeMixin&>(e).level = l;
		}
    };
	
public:
	
	/**
	 * This class works as a mixin by CRTP for instances of `Node` to be stored in 
	 * an spanning forest connected by the instances of `Edge`.
	 */
	template<class Node, class Edge>
	class HDTSpanningForestNodeMixin {
	
	public:
		
		HDTSpanningForestNodeMixin()
		: levelNodes() {
			levelNodes.push_back(new LevelNode<Node, Edge>(static_cast<Node&>(*this)));
		}
		
		HDTSpanningForestNodeMixin(const HDTSpanningForestNodeMixin<Node, Edge>&) = delete;
		
		HDTSpanningForestNodeMixin(HDTSpanningForestNodeMixin<Node, Edge>&&) = delete;
		
		~HDTSpanningForestNodeMixin() {
			for(LevelNode<Node, Edge>* n : levelNodes) {
				delete n;
			}
		}
		
	public:
		
		HDTSpanningForestNodeMixin<Node, Edge>& operator=(const HDTSpanningForestNodeMixin<Node, Edge>&) = delete;
		
		HDTSpanningForestNodeMixin<Node, Edge>& operator=(HDTSpanningForestNodeMixin<Node, Edge>&&) = delete;
		
	private:
		
		std::vector<LevelNode<Node, Edge>*> levelNodes;
		
		template<class, class> friend class HDTSpanningForestTrait;
	};
	
	/**
	 * This class works as a mixin by CRTP for instances of `Edge` to work as edges 
	 * in an spanning forest connecting the instances of `Node`.
	 */
	template<class Node, class Edge>
	class HDTSpanningForestEdgeMixin {
		
	public:
		
		HDTSpanningForestEdgeMixin()
		: node1(nullptr), node2(nullptr), level(0), isTreeEdge(false), levelEdges() {
		}
		
		HDTSpanningForestEdgeMixin(const HDTSpanningForestEdgeMixin<Node, Edge>&) = delete;
		
		HDTSpanningForestEdgeMixin(HDTSpanningForestEdgeMixin<Node, Edge>&&) = delete;
		
		~HDTSpanningForestEdgeMixin() {
			for(LevelEdge<Node, Edge>* e : levelEdges) {
				delete e;
			}
		}
		
	private:
		
		HDTSpanningForestEdgeMixin<Node, Edge>& operator=(const HDTSpanningForestEdgeMixin<Node, Edge>&) = delete;
		
		HDTSpanningForestEdgeMixin<Node, Edge>& operator=(HDTSpanningForestEdgeMixin<Node, Edge>&&) = delete;
		
	private:
		
		Node* node1;
		Node* node2;
		Level level;
		bool isTreeEdge;
		std::vector<LevelEdge<Node, Edge>*> levelEdges;
		
		template<class, class> friend class HDTSpanningForestTrait;
	};

	/**
	 * You can just use HDTSpanningForestNode and HDTSpanningForestEdge
	 * instead of writing classes extending HDTSpanningForestNodeMixin and HDTSpanningForestEdgeMixin,
	 * if you don't need to equip additional data or functionality.
	 */
	struct HDTSpanningForestNode;
	struct HDTSpanningForestEdge;

	struct HDTSpanningForestNode : public HDTSpanningForestNodeMixin<HDTSpanningForestNode, HDTSpanningForestEdge> {
	};

	struct HDTSpanningForestEdge : public HDTSpanningForestEdgeMixin<HDTSpanningForestNode, HDTSpanningForestEdge> {
	};
	
public:
	
	/*
	 * This provides operations on a spanning forest of `Node` and `Edge`.
	 * `Node` and `Edge` must implement the interface defined by the template class HDTSpanningForestTrait.
	 * Or you can provide an equivalent trait implementation by specifying the third template paramater.
	 */
	template<class Node, class Edge, class SFTrait = HDTSpanningForestTrait<Node, Edge>>
	class HDTSpanningForestAlgorithm {
		
	private:
		
		using LN = LevelNode<Node, Edge>;
		using LE = LevelEdge<Node, Edge>;
		using ETTAlgorithm = EulerTourTree::EulerTourTreeAlgorithm<LN, LE>;
		
	private:
	
		template<class CNode, class BackendItr>
		class ClusterIteratorImpl : public std::iterator<std::forward_iterator_tag, CNode> {
			
		private:
			
			using This = ClusterIteratorImpl<CNode, BackendItr>;
	
		public:
			
			ClusterIteratorImpl() 
			: itr() {
			}
			
			explicit ClusterIteratorImpl(BackendItr itr)
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
			
			CNode& operator*() const {
				return static_cast<CNode&>(itr->boxNode);
			}
			
			CNode* operator->() const {
				return static_cast<CNode*>(&itr->boxNode);
			}
		
		private:
			
			BackendItr itr;
		};
		
		template<class Trait>
		class EdgeIteratorImpl : public std::iterator<std::forward_iterator_tag, typename Trait::CEdge> {
			
		private:
			
			using This = EdgeIteratorImpl<Trait>;
			
		public:
			
			EdgeIteratorImpl()
			: node(nullptr), i(), j() {
			}
			
			EdgeIteratorImpl(typename Trait::CNode& n, typename Trait::CLevelNodeItr i)
			: node(&n), i(i), j() {
				if(i != SFTrait::levelNodes(*node).end()) {
					j = (*i)->edges.begin();
					stepForward();
				}
			}
			
		public:
			
			This& operator++() {
				++j;
				stepForward();
				return *this;
			}
			
			This operator++(int) {
				This tmp(*this);
				operator++();
				return tmp;
			}
			
			bool operator==(const This& other) const {
				return i == other.i && (i == SFTrait::levelNodes(*node).end() || j == other.j);
			}
			
			bool operator!=(const This& other) const {
				return i != other.i || (i != SFTrait::levelNodes(*node).end() && j != other.j);
			}
			
			typename Trait::CEdge& operator*() const {
				return **j;
			}
			
			typename Trait::CEdge* operator->() const {
				return *j;
			}
			
		private:
			void stepForward() {
				while(j == (*i)->edges.end()) {
					++i;
					if(i != SFTrait::levelNodes(*node).end()) {
						j = (*i)->edges.begin();
					} else {
						break;
					}
				}
			}
		
		private:
			
			typename Trait::CNode* node;
			typename Trait::CLevelNodeItr i;
			typename LN::EdgeSet::iterator j;
		};
		
	public:
		
		class Cluster {
			
		private:
			
			using Backend = typename ETTAlgorithm::NodeContainerView;
			
		private:
			
			explicit Cluster(Node& r)
			: rep(r), backend(ETTAlgorithm::nodeContainerView(*SFTrait::levelNodes(r)[0])) {
			}
			
		public:
			
			using iterator = ClusterIteratorImpl<Node, typename Backend::iterator>;
			using const_iterator = ClusterIteratorImpl<const Node, typename Backend::const_iterator>;
			
		public:
			
			iterator begin() {
				return iterator(backend.begin());
			}
			
			iterator end() {
				return iterator(backend.end());
			}
			
			const_iterator begin() const {
				return const_iterator(backend.cbegin());
			}
			
			const_iterator end() const {
				return const_iterator(backend.cend());
			}
			
			const_iterator cbegin() const {
				return const_iterator(backend.cbegin());
			}
			
			const_iterator cend() const {
				return const_iterator(backend.cend());
			}
			
			Node& representative() {
				return rep;
			}
			
			const Node& representative() const {
				return rep;
			}

			bool contains(const Node& n) const {
				return ETTAlgorithm::hasPath(*SFTrait::levelNodes(rep)[0], *SFTrait::levelNodes(n)[0]);
			}
			
			Size size() const {
				return ETTAlgorithm::clusterSize(*SFTrait::levelNodes(rep)[0]);
			}
		
		private:
			Node& rep;
			Backend backend;
			
			template<class, class, class> friend class HDTSpanningForestAlgorithm;
		};
		
		class ConstCluster {
			
		private:
			
			using Backend = typename ETTAlgorithm::NodeContainerView;
			
		private:
			
			explicit ConstCluster(const Node& r)
			: rep(r), backend(ETTAlgorithm::nodeContainerView(*SFTrait::levelNodes(r)[0])) {
			}
			
		public:
			
			using const_iterator = ClusterIteratorImpl<const Node, typename Backend::const_iterator>;
			
		public:
			
			const_iterator begin() const {
				return const_iterator(backend.cbegin());
			}
			
			const_iterator end() const {
				return const_iterator(backend.cend());
			}
			
			const_iterator cbegin() const {
				return const_iterator(backend.cbegin());
			}
			
			const_iterator cend() const {
				return const_iterator(backend.cend());
			}
			
			const Node& representative() const {
				return rep;
			}

			bool contains(const Node& n) const {
				return ETTAlgorithm::hasPath(*SFTrait::levelNodes(rep)[0], *SFTrait::levelNodes(n)[0]);
			}
			
			Size size() const {
				return ETTAlgorithm::clusterSize(*SFTrait::levelNodes(rep)[0]);
			}
		
		private:
			const Node& rep;
			Backend backend;
			
			template<class, class, class> friend class HDTSpanningForestAlgorithm;
		};
		
		class Edges {
			
		private:
			
			explicit Edges(Node& n)
			: node(n) {
			}
			
		private:
			
			struct T {
				using CNode = Node;
				using CEdge = Edge;
				using CLevelNodeItr = typename std::vector<LN*>::iterator;
			};
			
			struct CT {
				using CNode = const Node;
				using CEdge = const Edge;
				using CLevelNodeItr = typename std::vector<LN*>::const_iterator;
			};
			
		public:
			
			using iterator = EdgeIteratorImpl<T>;
			using const_iterator = EdgeIteratorImpl<CT>;
			
		public:
			
			iterator begin() {
				return iterator(node, SFTrait::levelNodes(node).begin());
			}
			
			iterator end() {
				return iterator(node, SFTrait::levelNodes(node).end());
			}
			
			const_iterator begin() const {
				return const_iterator(node, SFTrait::levelNodes(node).begin());
			}
			
			const_iterator end() const {
				return const_iterator(node, SFTrait::levelNodes(node).end());
			}
			
			const_iterator cbegin() const {
				return const_iterator(node, SFTrait::levelNodes(node).begin());
			}
			
			const_iterator cend() const {
				return const_iterator(node, SFTrait::levelNodes(node).end());
			}
			
		private:
			Node& node;
			
			template<class, class, class> friend class HDTSpanningForestAlgorithm;
		};
		
		class ConstEdges {
			
		private:
			
			explicit ConstEdges(const Node& n)
			: node(n) {
			}
			
		private:
			
			struct CT {
				using CNode = const Node;
				using CEdge = const Edge;
				using CLevelNodeItr = typename std::vector<LN*>::const_iterator;
			};
			
		public:
			
			using const_iterator = EdgeIteratorImpl<CT>;
			
			
		public:
			
			const_iterator begin() const {
				return const_iterator(node, SFTrait::levelNodes(node).cbegin());
			}
			
			const_iterator end() const {
				return const_iterator(node, SFTrait::levelNodes(node).cend());
			}
			
			const_iterator cbegin() const {
				return const_iterator(node, SFTrait::levelNodes(node).cbegin());
			}
			
			const_iterator cend() const {
				return const_iterator(node, SFTrait::levelNodes(node).cend());
			}
			
		private:
			const Node& node;
			
			template<class, class, class> friend class HDTSpanningForestAlgorithm;
		};
		
		
	public:
		
		/**
		 * Returns whether or not `n1` and `n2` are in the same connected component.
		 */
		static bool hasPath(const Node& n1, const Node& n2) {
			return ETTAlgorithm::hasPath(
				*SFTrait::levelNodes(n1)[0],
				*SFTrait::levelNodes(n2)[0]
			);
		}
		
		/**
		 * Returns a Cluster which represents the connected compnent of `n`.
		 * It has an interface for STL compatible iterators with which
		 * you can visit each node in the component one by one.
		 */
		static Cluster cluster(Node& n) {
			return Cluster(
				ETTAlgorithm::findClusterRep(
					*SFTrait::levelNodes(n)[0]
				).boxNode
			); 
		}
		
		/**
		 * Returns a ConstCluster which represents the connected compnent of `n`.
		 * It has an interface for STL compatible iterators with which
		 * you can visit each node in the component one by one.
		 */
		static ConstCluster cluster(const Node& n) {
			return ConstCluster(
				ETTAlgorithm::findClusterRep(
					*SFTrait::levelNodes(n)[0]
				).boxNode
			); 
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
			return ETTAlgorithm::clusterSize(
				ETTAlgorithm::findClusterRep(
					*SFTrait::levelNodes(n)[0]
				)
			);
		}
		
		/**
		 * Returns whether of not `n` is the representative of its connected component.
		 * Each connected component has a representative.
		 * The representative can be changed by any operations modifying the connected component.
		 */
		static bool isClusterRep(const Node& n) {
			return &ETTAlgorithm::findClusterRep(
				*SFTrait::levelNodes(n)[0]
			).boxNode == &n;
		}
		
		/**
		 * Returns the representative of the connected component of `n`.
		 */
		static const Node& findClusterRep(const Node& n) {
			return static_cast<const Node&>(
				ETTAlgorithm::findClusterRep(
					*SFTrait::levelNodes(n)[0]
				).boxNode
			);
		}
		
		/**
		 * Returns the representative of the connected component of `n`.
		 */
		static Node& findClusterRep(Node& n) {
			return static_cast<Node&>(
				ETTAlgorithm::findClusterRep(
					*SFTrait::levelNodes(n)[0]
				).boxNode
			);
		}
		
		/**
		 * Returns an Edges with which you can visit the edges of `n` one by one.
		 * It has an interface for STL compatible iterators. 
		 */
		static Edges edges(Node& n) {
			return Edges(n);
		}
		
		/**
		 * Returns an ConstEdges with which you can visit the edges of `n` one by one.
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
			SFTrait::level(e, 0);
			SFTrait::levelNodes(n1)[0]->edges.insert(&e);
			SFTrait::levelNodes(n2)[0]->edges.insert(&e);
			if(ETTAlgorithm::hasPath(
				*SFTrait::levelNodes(n1)[0], 
				*SFTrait::levelNodes(n2)[0]
			)) {
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
			if(SFTrait::isTreeEdge(e)) {
				eraseTreeEdge(e);
				clusterSplitted = ! checkReplacement(e);
			} else {
				eraseNonTreeEdge(e);
			}
			SFTrait::invalidate(e);
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
		
		static void eraseTreeEdge(Edge& e) {
			for(LE* ee : SFTrait::levelEdges(e)) {
				ETTAlgorithm::disconnect(*ee);
				delete ee;
			}
			SFTrait::levelEdges(e).clear();
			levelNode1(e, SFTrait::level(e)).edges.erase(&e);
			levelNode2(e, SFTrait::level(e)).edges.erase(&e);
		}
		
		static void eraseNonTreeEdge(Edge& e) {
			levelNode1(e, SFTrait::level(e)).edges.erase(&e);
			levelNode2(e, SFTrait::level(e)).edges.erase(&e);
		}
		
		static bool checkReplacement(Edge& e) {
			for(Level l = 0; l <= SFTrait::level(e); ++l) {
				if(checkReplacement(e, SFTrait::level(e) - l)) {
					return true;
				}
			}
			return false;
		}
		
		static bool checkReplacement(Edge& e, Level l) {
			LN& r1 = ETTAlgorithm::findClusterRep(levelNode1(e, l));
			LN& r2 = ETTAlgorithm::findClusterRep(levelNode2(e, l));
			typename ETTAlgorithm::Size sz1 = ETTAlgorithm::clusterSize(r1);
			typename ETTAlgorithm::Size sz2 = ETTAlgorithm::clusterSize(r2);
			if(sz1 < sz2) {
				return checkReplacement(r1, r2, l);
			} else {
				return checkReplacement(r2, r1, l);
			}
		}

		static bool checkReplacement(LN& smaller, LN& larger, Level l) {
			levelUpTreeEdges(smaller, l);
			LN& largerRep = ETTAlgorithm::findClusterRep(larger);
			for(LN& n : ETTAlgorithm::nodeContainerView(smaller)) {
				auto i = n.edges.begin(); 
				while(i != n.edges.end()) {
					Edge& e = **i;
					LN& m = &levelNode1(e, l) == &n ? levelNode2(e, l) : levelNode1(e, l);
					if(&ETTAlgorithm::findClusterRep(m) == &largerRep) {
						replaceWith(e);
						return true;
					} else {
						i = n.edges.erase(i);
						m.edges.erase(&e);
						putInUpperLevel(e, l);
					}
				}
			}
			return false;
		}

		static void levelUpTreeEdges(LN& n, Level l) {
			for(LE& le : ETTAlgorithm::edgeContainerView(n)) {
				Edge& e = le.boxEdge;
				if(SFTrait::level(e) == l) {
					//level up the box edge
					levelNode1(e, l).edges.erase(&e);
					levelNode2(e, l).edges.erase(&e);
					putInUpperLevel(e, l);
					//add one more level edge
					LE* newLE = new LE(e);
					SFTrait::levelEdges(e).push_back(newLE);
					ETTAlgorithm::connect(levelNode1(e, l + 1), levelNode2(e, l + 1), *newLE);
				}
			}
		}

		static void putInUpperLevel(Edge& e, Level l) {
			SFTrait::level(e, l + 1);
			putInUpperLevel(e, SFTrait::node1(e), l);
			putInUpperLevel(e, SFTrait::node2(e), l);
		}

		static void putInUpperLevel(Edge& e, Node& n, Level l) {
			std::vector<LevelNode<Node, Edge>*>& vec = SFTrait::levelNodes(n);
			if(vec.size() < l + 2) {
				vec.push_back(new LN(n));
			}
			vec[l + 1]->edges.insert(&e);
		}

		static void replaceWith(Edge& e) {
			SFTrait::levelEdges(e).reserve(SFTrait::level(e) + 1);
			for(Level l = 0; l <= SFTrait::level(e); ++l) {
				LE* le = new LE(e);
				SFTrait::levelEdges(e).push_back(le);
				ETTAlgorithm::connect(levelNode1(e, l), levelNode2(e, l), *le);
			}
			SFTrait::setTreeEdge(e, true);
		}

		static LN& levelNode1(Edge& e, Level l) {
			return *SFTrait::levelNodes(SFTrait::node1(e))[l];
		}

		static LN& levelNode2(Edge& e, Level l) {
			return *SFTrait::levelNodes(SFTrait::node2(e))[l];
		}
	};
	
};


/**
 * Public interface from HDTSpanningForestModule.
 */

template<class Node, class Edge>
using HDTSpanningForestTrait = HDTSpanningForestModule::HDTSpanningForestTrait<Node, Edge>;

template<class Node, class Edge>
using HDTSpanningForestAlgorithm = HDTSpanningForestModule::HDTSpanningForestAlgorithm<Node, Edge>;

template<class Node, class Edge>
using HDTSpanningForestNodeMixin = HDTSpanningForestModule::HDTSpanningForestNodeMixin<Node, Edge>;

template<class Node, class Edge>
using HDTSpanningForestEdgeMixin = HDTSpanningForestModule::HDTSpanningForestEdgeMixin<Node, Edge>;

using HDTSpanningForestNode = HDTSpanningForestModule::HDTSpanningForestNode;

using HDTSpanningForestEdge = HDTSpanningForestModule::HDTSpanningForestEdge;

} //end namespace HDTSpanningForest
} //end namespace Cnrc
} //end namespace Snu
