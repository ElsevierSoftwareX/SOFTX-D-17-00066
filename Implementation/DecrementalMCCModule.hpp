/**
 * An implementation of the decremental dynamic algorithm for mutually connected components 
 * by S. Hwang, S. Choi, D. Lee and B. Kahng (Physical Review E, 91, 022814 (2015)).
 */

#pragma once

#include <queue>
#include <set>
#include <map>
#include <vector>
#include <random>
#include "EulerTourTreeSpanningForest.hpp"

namespace Snu {
namespace Cnrc {
namespace DecrementalMCCModule {

template<
	template<class, class> class SpanningForestNodeMixin,
	template<class, class> class SpanningForestEdgeMixin,
	template<class, class> class SpanningForestAlgorithm
>
class DecrementalMCCModule {

public:

	template <class, class> class DecrementalMCCNodeMixin;
	template <class, class> class DecrementalMCCEdgeMixin;
	template <class, class> class DecrementalMCCAlgorithm;

	using Size = unsigned int;

private:

	using NetworkId = int;	
	static const NetworkId DEFAULT_NETID = -1;
	static const NetworkId FIRST_NETID = 0;
	static const NetworkId SECOND_NETID = 1;

public:

	/**
	 * This trait defines an interface for the type `Node` and `Edge` to be manipulated by the operations
	 * provided by DecrementalMCCAlgorithm<Node, Edge>.
	 * If `Node` and `Edge` extend DecrementalMCCNodeMixin<Node, Edge> and DecrementalMCCEdgeMixin<Node, Edge> respectively,
	 * then the default implementation will just work.
	 * If not, a specialization of this template for `Node` and `Edge` must be provided.
	 */
	template<class Node, class Edge>
	class DecrementalMCCTrait {

	private:

		using NodeMixin = DecrementalMCCNodeMixin<Node, Edge>;
		using EdgeMixin = DecrementalMCCEdgeMixin<Node, Edge>;

		static_assert(
			std::is_base_of<NodeMixin, Node>::value, 
			"The default implementation of DecrementalMCCTrait<S, T> requires "
			"S extending DecrementalMCCNodeMixin<S, T>. "
			"If you don't want this, you need to provide a specialization of DecrementalMCCTrait for S, T."
		);

		static_assert(
			std::is_base_of<EdgeMixin, Edge>::value, 
			"The default implementation of DecrementalMCCTrait<S, T> requires "
			"T extending DecrementalMCCEdgeMixin<S, T>. "
			"If you don't want this, you need to provide a specialization of DecrementalMCCTrait for S, T."
		);

	public:

		static std::vector<Node*>& interNeighbors(Node& n) {
			return static_cast<NodeMixin&>(n).interNeighbors;
		}

		static const std::vector<Node*>& interNeighbors(const Node& n) {
			return static_cast<const NodeMixin&>(n).interNeighbors;
		}

		static int networkId(const Node& n) {
			return static_cast<const NodeMixin&>(n).networkId;
		}

		static void networkId(Node& n, int id) {
			static_cast<NodeMixin&>(n).networkId = id;
		}

		static bool isAlienated(const Node& n) {
			return static_cast<const NodeMixin&>(n).toAlienate;
		}

		static void setAlienated(Node& n, bool b) {
			static_cast<NodeMixin&>(n).toAlienate = b;
		}
	};

	/**
	 * This class works as a mixin by CRTP for instances of `Node` to be manipulated by
	 * the operations of DecrementalMCCAlgorithm<Node, Edge>.
	 */
	template<class Node, class Edge>
	class DecrementalMCCNodeMixin : public SpanningForestNodeMixin<Node, Edge> {

	public:
		DecrementalMCCNodeMixin()
		: networkId(DEFAULT_NETID), toAlienate(false), interNeighbors() {
		}

	private:
		NetworkId networkId;
		bool toAlienate;
		std::vector<Node*> interNeighbors;

		template<class, class> friend class DecrementalMCCTrait;
	};

	/**
	 * This class works as a mixin by CRTP for instances of `Edge` to be manipulated by
	 * the operations of DecrementalMCCAlgorithm<Node, Edge>.
	 */
	template<class Node, class Edge>
	class DecrementalMCCEdgeMixin : public SpanningForestEdgeMixin<Node, Edge> {
	};

	
	/**
	 * You can just use DecrementalMCCNode and DecrementalMCCEdge
	 * instead of writing classes extending DecrementalMCCNodeMixin and DecrementalMCCEdgeMixin,
	 * if you don't need to equip additional data or functionality.
	 */
	class DecrementalMCCNode;
	class DecrementalMCCEdge;

	class DecrementalMCCNode : public DecrementalMCCNodeMixin<DecrementalMCCNode, DecrementalMCCEdge> {
	};

	class DecrementalMCCEdge : public DecrementalMCCEdgeMixin<DecrementalMCCNode, DecrementalMCCEdge> {
	};


	/**
	 * This provides operations to trace the MCCs in interdependent networks whose node and edges are 
	 * represented by `Node` and `Edge` respectively during the nodes or edges are deleted one by one.
	 */	
	//template<class Node, class Edge, class Trait = DecrementalMCCTrait<Node, Edge>>
	template<class Node, class Edge>
	class DecrementalMCCAlgorithm {

	private:
		/* 
		 * IMPLEMENTATION NOTE:
		 * Putting this as the third template parameter with default is more flexible.
		 * However, GCC 6.4/7.1 fails to compile it while clang 4.0 and VC++ 2017 successfully compile it.
		 * This workaround is less flexible than the template parameter way in strict sense
		 * but makes no difference in practice.
		 */
		using Trait = DecrementalMCCTrait<Node, Edge>;

	public:
		using SpanningForest = SpanningForestAlgorithm<Node, Edge>;

	private:
		struct SortByClusterSize {
			bool operator()(const Node* n, const Node* m) const {
				const Size sz1 = SpanningForest::clusterSize(*n);
				const Size sz2 = SpanningForest::clusterSize(*m);
				return sz1 == sz2 ? n > m : sz1 > sz2;
			}
		};

		struct PhantomClusterRepresentativeSet {
			void insert(Node*) {
			}
			
			std::size_t erase(Node*) {
				return 0;
			}
		};

		struct PhantomClusterSizeDist {
			Size& operator[](Size) {
				return p.second;
			}
			
			void erase(std::pair<Size, Size>*) {
			}

			std::pair<Size, Size>* begin() {
				return &p;
			}

			std::pair<Size, Size>* end() {
				return &p;
			}

			std::pair<Size, Size>* find(Size) {
				return &p;
			}

			std::pair<Size, Size> p;
			Size alienatedCount = 0;
		};

	public:
		class ClusterSizeDist : public std::map<Size, Size> {

		public:

			using std::map<Size, Size>::map;

		private:

			Size alienatedCount = 0;

			template<class, class> friend class DecrementalMCCAlgorithm;
		};

		class ClusterRepresentativeSet : public std::set<const Node*, SortByClusterSize> {
		};

	public:
		/**
		 * Set an interconnection between `n1` and `n2`.
		 * The nodes `n1` and `n2` should be in different networks.
		 */
		static void interconnect(Node& n1, Node& n2) {
			Trait::interNeighbors(n1).push_back(&n2);
			Trait::interNeighbors(n2).push_back(&n1);
		}

		/**
		 * Returns the `std::vector` of `Node` pointers interconnected with `n`.
		 */
		static const std::vector<Node*>& interconnectedNeighbors(const Node& n) {
			return Trait::interNeighbors(n);
		}

		/**
		 * Disconnects the edge `e` and updates `cdist`.
		 * See the comment on `initialize` for the details.
		 */
		static void disconnect(Edge& e, ClusterSizeDist& cdist1) {
			PhantomClusterRepresentativeSet reprs;
			disconnect_(e, cdist1, reprs, reprs);
		}

		/**
		 * Disconnects the edge `e` and updates `cdist`, `reprs1`.
		 * See the comment on `initialize` for the details.
		 */
		static void disconnect(Edge& e, ClusterSizeDist& cdist1, ClusterRepresentativeSet& reprs1) {
			PhantomClusterRepresentativeSet reprs2;
			disconnect_(e, cdist1, reprs1, reprs2);
		}

		/**
		 * Disconnects the edge `e` and updates `reprs1` and `reprs2`.
		 * See the comment on `initialize` for the details.
		 */
		static void disconnect(Edge& e, ClusterRepresentativeSet& reprs1, ClusterRepresentativeSet& reprs2) {
			PhantomClusterSizeDist cdist1;
			disconnect_(e, cdist1, reprs1, reprs2);
		}

		/**
		 * Isolates the node `n` from its network.
		 * See the comment on `initialize` for the details.
		 */
		static void disconnect(Node& n, ClusterSizeDist& cdist1) {
			PhantomClusterRepresentativeSet reprs;
			disconnect_(n, cdist1, reprs, reprs);
		}

		/**
		 * Isolates the node `n` from its network
		 * See the comment on `initialize` for the details.
		 */
		static void disconnect(Node& n, ClusterSizeDist& cdist1, ClusterRepresentativeSet& reprs1) {
			ClusterRepresentativeSet reprs2;
			disconnect_(n, cdist1, reprs1, reprs2);
		}

		/**
		 * Isolates the node `n` from its network
		 * See the comment on `initialize` for the details.
		 */
		static void disconnect(Node& n, ClusterRepresentativeSet& reprs1, ClusterRepresentativeSet& reprs2) {
			PhantomClusterSizeDist cdist1;
			disconnect_(n, cdist1, reprs1, reprs2);
		}

		/**
		 * NodeContainer must contains instances of `Node` and be iterable by the range-based for loop and also provide `size` method.
		 * ClusterSizeDist has the same interface with std::map<Size, Size> and `cdist` must be empty.
		 *
		 * All `Node` instances in `nodes1` will be assumed to be in the same network in the subsequent calls of `disconnect` method.
		 * We call the network of the nodes in `nodes1` the first network.
		 * The nodes in `nodes2` are also assumed to be in a single network and we call it the second network.
		 * You can create edges in each network using `SpanningForest::connect` method and set the interdependency between the two networks
		 * using the `interconnect` method. Those configurations must be done before calling this method.
		 *
		 * For each MCC, there is a corresponding connected component in each network.
		 * We call the connected component in a network as the projected MCC on the network.
		 *
		 * The call of `initialize` updates `cdist` so that it contains the size distribution of MCCs **projectced on the first network**.
		 * If each node has one and only one interconnection, the second network also has the same distribution.
		 * Otherwise, however, this is not guaranteed. 
		 * Even some nodes may not belong to any MCC and thus the projected components don't partition the network.
		 * In this case, you need to use another version of `initialize` although it costs more space and time.
		 *
		 * To trace the size distribution dynamically after initialization by this method, call any version of `disconnect` method
		 * with the same `cdist` given to this method.
		 */
		template<class NodeContainer>
		static void initialize(NodeContainer& nodes1, NodeContainer& nodes2, ClusterSizeDist& cdist1) {
			PhantomClusterRepresentativeSet reprs;
			initialize_(nodes1, nodes2, cdist1, reprs, reprs);
		}

		/**
		 * See the comment on the other version of `initialize` first.
		 *
		 * This method stores the representative nodes of MCCs projected into the first network at `reprs1`.
		 * You can query the size of the component and iterate over the nodes in the component, and so on with the representatives.
		 * For example, `SpanningForest::findClusterRep(interconnectedNeighbors(n)[0])` returns the representative of the component interconnected to `n` 
		 * (i.e., the representative of the connected component which is the projection of the MCC of the node `n` into the other network).
		 *
		 * ClusterRepresentativeSet has the same interface with std::set<Node*> and is sorted by the size of the component in decreasing order.
		 * Thus, e.g., `*reprs1.begin()` is the pointer to the representative of the largest projected MCC in the first network.
		 * You must give empty `reprs1` to this method.
		 *
		 * To trace the size distribution dynamically after initialization by this method, call any version of `disconnect` method
		 * with the same `cdist` and `reprs1` given to this method.
		 */
		template<class NodeContainer>
		static void initialize(
			NodeContainer& nodes1, 
			NodeContainer& nodes2, 
			ClusterSizeDist& cdist1,
			ClusterRepresentativeSet& reprs1
		) {
			PhantomClusterRepresentativeSet reprs2;
			initialize_(nodes1, nodes2, cdist1, reprs1, reprs2);
		}

		/**
		 * See the comment on the other two versions of `initialize` first.
		 * This method stores the representative nodes of MCCs projected on the first networks at `reprs1` and on the second network at `reprs2`.
		 */
		template<class NodeContainer>
		static void initialize(
			NodeContainer& nodes1, 
			NodeContainer& nodes2, 
			ClusterRepresentativeSet& reprs1,
			ClusterRepresentativeSet& reprs2
		) {
			PhantomClusterSizeDist cdist1;
			initialize_(nodes1, nodes2, cdist1, reprs1, reprs2);
		}

	private:

		template<class NodeContainer, class CDist, class ReprSet1, class ReprSet2>
		static void initialize_(
			NodeContainer& nodes1, 
			NodeContainer& nodes2, 
			CDist& cdist,
			ReprSet1& reprs1,
			ReprSet2& reprs2
		) {
			for(Node& n : nodes1) {
				Trait::networkId(n, FIRST_NETID);
				Trait::setAlienated(n, false);
			}
			for(Node& n : nodes2) {
				Trait::networkId(n, SECOND_NETID);
				Trait::setAlienated(n, false);
			}
			std::vector<Edge> adhocEdges1 = createAdhocEdges(nodes1);
			std::vector<Edge> adhocEdges2 = createAdhocEdges(nodes2);
			cdist[nodes1.size()] = 1;
			reprs1.insert(&SpanningForest::findClusterRep(*nodes1.begin()));
			reprs2.insert(&SpanningForest::findClusterRep(*nodes2.begin()));
			for(Node& n : nodes1) {
				if(Trait::interNeighbors(n).size() == 0) {
					disconnect_(n, cdist, reprs1, reprs2);
				}
			}
			for(Node& n : nodes2) {
				if(Trait::interNeighbors(n).size() == 0) {
					disconnect_(n, cdist, reprs1, reprs2);
				}
			}
			for(Edge& e : adhocEdges1) {
				disconnect_(e, cdist, reprs1, reprs2);
			}
			for(Edge& e : adhocEdges2) {
				disconnect_(e, cdist, reprs1, reprs2);
			}
		}

		template<class NodeContainer>
		static std::vector<Edge> createAdhocEdges(NodeContainer& nodes) {
			std::set<Node*> reprs;
			for(Node& n : nodes) {
				reprs.insert(&SpanningForest::findClusterRep(n));
			}
			std::vector<Edge> adhocEdges(reprs.size() - 1);
			auto i2 = reprs.begin();
			auto i1 = i2++;
			Size i0 = 0;
			while(i2 != reprs.end()) {
				SpanningForest::connect(**i1, **i2, adhocEdges[i0]);
				++i0;
				++i1;
				++i2;
			}
			return adhocEdges;//RVO or move
		}

		template<class CDist>
		static void adjustClusterSizeDist(CDist& cdist) {
			auto i = cdist.find(1);
			if(i != cdist.end()) {
				i->second -= cdist.alienatedCount;
				if(i->second == 0) {
					cdist.erase(i);
				}
			}
		}

		template<class CDist>
		static void putClusterSizes(CDist& cdist, Size sz1, Size sz2) {
			cdist[sz1]++;
			cdist[sz2]++;
			auto i = cdist.find(sz1 + sz2);
			if(i->second == 1) {
				cdist.erase(i);
			} else {
				i->second--;
			}
		}

		template<class CDist, class ReprSet1, class ReprSet2>
		static void disconnect_(Edge& e, CDist& cdist, ReprSet1& reprs1, ReprSet2& reprs2) {
			std::queue<Edge*> toInactivate;
			toInactivate.push(&e);
			std::queue<Node*> toAlienate;
			cdist.alienatedCount = 0;
			disconnect_(toInactivate, toAlienate, cdist, reprs1, reprs2);
			adjustClusterSizeDist(cdist);
		}

		template<class CDist, class ReprSet1, class ReprSet2>
		static void disconnect_(Node& n, CDist& cdist, ReprSet1& reprs1, ReprSet2& reprs2) {
			std::queue<Edge*> toInactivate;
			for(Edge& e : SpanningForest::edges(n)) {
				toInactivate.push(&e);
			}
			std::queue<Node*> toAlienate;
			if(Trait::interNeighbors(n).size() == 0) {
				toAlienate.push(&n);
			}
			cdist.alienatedCount = 0;
			disconnect_(toInactivate, toAlienate, cdist, reprs1, reprs2);
			adjustClusterSizeDist(cdist);
		}

		template<class CDist, class ReprSet1, class ReprSet2>
		static void disconnect_(
			std::queue<Edge*>& toInactivate,
			std::queue<Node*>& toAlienate,
			CDist& cdist,
			ReprSet1& reprs1,
			ReprSet2& reprs2
		) {
			while(! toInactivate.empty() || ! toAlienate.empty()) {
				if(! toInactivate.empty()) {
					Edge* ti = toInactivate.front();
					toInactivate.pop();
					if(SpanningForest::isValid(*ti)) {
						if(Trait::networkId(SpanningForest::node1(*ti)) == FIRST_NETID) {
							disconnect_(*ti, toInactivate, toAlienate, cdist, reprs1);
						} else {
							disconnect_(*ti, toInactivate, toAlienate, cdist, reprs2);
						}
					}
				}
				if(! toAlienate.empty()) {
					Node* an = toAlienate.front();
					toAlienate.pop();
					if(Trait::networkId(*an) == FIRST_NETID) {
						alienate(*an, toInactivate, toAlienate, cdist, reprs1);
					} else {
						alienate(*an, toInactivate, toAlienate, cdist, reprs2);
					}
				}
			}
		}

		template<class CDist, class ReprSet>
		static void disconnect_(
			Edge& ti,
			std::queue<Edge*>& toInactivate,
			std::queue<Node*>& toAlienate,
			CDist& cdist,
			ReprSet& reprs
		) {
			Node& n1 = SpanningForest::node1(ti);
			Node& n2 = SpanningForest::node2(ti);
			Node& oldRoot = SpanningForest::findClusterRep(n1);
			reprs.erase(&oldRoot);//must be called before the disconnect call
			if(! SpanningForest::disconnect(ti)) {
				Node& newRoot = SpanningForest::findClusterRep(n1);
				reprs.insert(&newRoot);
				return;
			}
			using Clst = typename SpanningForest::Cluster;
			Clst c1 = SpanningForest::cluster(n1);
			Clst c2 = SpanningForest::cluster(n2);
			Clst smaller = c1.size() < c2.size() ? c1 : c2;
			Clst larger  = c1.size() < c2.size() ? c2 : c1;
			if(FIRST_NETID == Trait::networkId(n1)) {
				putClusterSizes(cdist, c1.size(), c2.size());
			}
			for(Node& nInA : smaller) {
				checkInconsistency(nInA, larger.representative(), toInactivate, toAlienate);
			}
			if(c1.size() != 1 || ! Trait::isAlienated(n1)) {
				reprs.insert(&c1.representative());
			}
			if(c2.size() != 1 || ! Trait::isAlienated(n2)) {
				reprs.insert(&c2.representative());
			}
		}

		static void checkInconsistency(
			Node& nInA,
			Node& largerClusterRepr,
			std::queue<Edge*>& toInactivate,
			std::queue<Node*>& toAlienate
		) {
			for(Node* nInB : Trait::interNeighbors(nInA)) {
				for(Node* xnInA : Trait::interNeighbors(*nInB)) {
					if(&SpanningForest::findClusterRep(*xnInA) == &largerClusterRepr) {
						toAlienate.push(nInB);
						break;
					}
				}
				for(Edge& eInB : SpanningForest::edges(*nInB)) {
					Node& mInB = &SpanningForest::node1(eInB) == nInB ? 
								SpanningForest::node2(eInB): 
								SpanningForest::node1(eInB);
					for(Node* mInA : Trait::interNeighbors(mInB)) {
						if(&SpanningForest::findClusterRep(*mInA) == &largerClusterRepr) {
							toInactivate.push(&eInB);
						}
					}
				}
			}
		}

		template<class CDist, class ReprSet>
		static void alienate(
			Node& an,
			std::queue<Edge*>& toInactivate,
			std::queue<Node*>& toAlienate,
			CDist& cdist,
			ReprSet& reprs
		) {
			if(Trait::isAlienated(an)) {
				return;
			}
			Trait::setAlienated(an, true);
			if(Trait::networkId(an) == FIRST_NETID) {
				++cdist.alienatedCount;
			}
			if(SpanningForest::edges(an).begin() != SpanningForest::edges(an).end()) {
				for(Edge& e : SpanningForest::edges(an)) {
					toInactivate.push(&e);
				}
			} else {
				reprs.erase(&an); //if any
			}
			for(Node* nInB : Trait::interNeighbors(an)) {
				toAlienate.push(nInB);
			}
		}

	};

};

} // end namespace DecrementalMCCModule
} // end namespace Cnrc
} // end namespace Snu
