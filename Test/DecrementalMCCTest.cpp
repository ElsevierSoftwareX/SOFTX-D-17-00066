#define CATCH_CONFIG_MAIN
#include "catch/catch.hpp"
#include <vector>
#include <map>
#include <set>
#include <utility>

#ifdef CNRC_ENABLE_HDT
	#include "DecrementalMCCWithHDT.hpp"
	namespace MCC = Snu::Cnrc::DecrementalMCCWithHDT;
#else
	#include "DecrementalMCC.hpp"
	namespace MCC = Snu::Cnrc::DecrementalMCC;
#endif

using namespace std;

struct Node;
struct Edge;
struct Node : MCC::DecrementalMCCNodeMixin<Node, Edge> {
	template<class T>
	T& operator<<(T& out) {
		out << name;
	}
	int name;
};
struct Edge : MCC::DecrementalMCCEdgeMixin<Node, Edge> {
};

using NodeVector = vector<Node>;
using EdgeVector = vector<Edge>;
using DecrementalMCCAlgorithm = MCC::DecrementalMCCAlgorithm<Node, Edge>;
using SpanningForest = DecrementalMCCAlgorithm::SpanningForest;
using ClusterSizeDist = DecrementalMCCAlgorithm::ClusterSizeDist;
using ClusterRepresentativeSet = DecrementalMCCAlgorithm::ClusterRepresentativeSet;
using Size = MCC::Size;

set<const Node*> cluster(const Node& n) {
	set<const Node*> s;
	for(const Node& m : SpanningForest::cluster(n)) {
		REQUIRE(s.insert(&m).second);
	}
	return s;
}

template<class...Args>
set<const Node*> nodeSet(const NodeVector& nodes, Args...args) {
	vector<int> sv{args...};
	set<const Node*> ret;
	for(auto i : sv) {
		ret.insert(&nodes[i]);
	}
	return ret;
}

bool checkReprSet(const NodeVector& nodes, const ClusterRepresentativeSet& reprs, size_t sz) {
	using Trait = MCC::DecrementalMCCTrait<Node, Edge>;
	set<const Node*> checkedReprs;
	ClusterRepresentativeSet mccReprs;
	size_t largestSize = 0;
	for(const Node& n : nodes) {
		SpanningForest::ConstCluster c = SpanningForest::cluster(n);
		if(checkedReprs.count(&c.representative())) {
			continue;
		}
		checkedReprs.insert(&c.representative());
		bool b = true;
		for(const Node& m : c) {
			b &= Trait::interNeighbors(m).size() > 0;
		}
		if(! b) {
			continue;
		}
		SpanningForest::ConstCluster d = SpanningForest::cluster(
			static_cast<const Node&>(*Trait::interNeighbors(n)[0])
		);
		for(const Node& m : d) {
			b &= Trait::interNeighbors(m).size() > 0;
		}
		for(const Node& m : c) {
			for(const Node* l : Trait::interNeighbors(m)) {
				b &= d.contains(*l);
			}
		}
		for(const Node& m : d) {
			for(const Node* l : Trait::interNeighbors(m)) {
				b &= c.contains(*l);
			}
		}
		if(b) {
			mccReprs.insert(&c.representative());
		}
		if(largestSize < c.size()) {
			largestSize = c.size();
		}
	}
	return reprs.size() == sz && mccReprs == reprs && 
		(sz == 0 || largestSize == SpanningForest::clusterSize(**reprs.begin()));
}

bool ensureTotallyDisconnected(const NodeVector& nodes) {
	bool b = true;
	for(size_t i = 0; i < nodes.size(); ++i) {
		b &= cluster(nodes[i]) == nodeSet(nodes, static_cast<int>(i));
	}
	for(size_t i = 0; i < nodes.size(); ++i) {
		for(size_t j = 0; j < nodes.size(); ++j) {
			if(i != j) {
				b &= ! SpanningForest::hasPath(nodes[i], nodes[j]);
			}
		}
	}
	return b;
}

bool ensureTotallyDisconnected(const NodeVector& nodes1, const NodeVector& nodes2) {
	return ensureTotallyDisconnected(nodes1) && ensureTotallyDisconnected(nodes2);
}

TEST_CASE("TestSingleNodeInterdependentNetworks-SingleInterdepedency") {
	// Trivial single node networks
	NodeVector nodes1(1);
	NodeVector nodes2(1);
	EdgeVector edges1;
	EdgeVector edges2;

	// Set interconnections between the two networks.
	DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[0]);
	REQUIRE(DecrementalMCCAlgorithm::interconnectedNeighbors(nodes1[0]).front() == nodes2.data());

	SECTION("") {
		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));
	}

	SECTION("") {
		// Initial distribution of MCC sizes.
		ClusterSizeDist cdist1;
		ClusterRepresentativeSet reprs1;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, cdist1, reprs1);
		REQUIRE(cdist1 == (ClusterSizeDist{{1, 1}}));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
	}

	SECTION("") {
		// Initial distribution of MCC sizes.
		ClusterSizeDist cdist1;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, cdist1);
		REQUIRE(cdist1 == (ClusterSizeDist{{1, 1}}));
	}
}

TEST_CASE("TestTwoNodesInterdependentNetworks-SingleInterdepedency") {

	NodeVector nodes1(2); // nodes of network #1
	NodeVector nodes2(2); // nodes of network #2
	EdgeVector edges1(1); // edges of network #1
	EdgeVector edges2(1); // edges of network #2

	// Set interconnections between the two networks.
	// Each node has one interconnection.
	DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[0]);
	DecrementalMCCAlgorithm::interconnect(nodes1[1], nodes2[1]);
	
	SECTION("") {
		// Each network contains no edges at all.

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(checkReprSet(nodes1, reprs1, 2));
		REQUIRE(checkReprSet(nodes2, reprs2, 2));
	}

	SECTION("") {
		// Each network contains no edges at all.

		// Initial distribution of MCC sizes.
		ClusterSizeDist cdist1;
		ClusterRepresentativeSet reprs1;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, cdist1, reprs1);
		REQUIRE(cdist1 == (ClusterSizeDist{{1, 2}}));
		REQUIRE(checkReprSet(nodes1, reprs1, 2));
	}

	SECTION("") {
		// Configuration of network #1. Network #2 has no edges.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(! SpanningForest::hasPath(nodes1[0], nodes1[1]));
		REQUIRE(! SpanningForest::hasPath(nodes2[0], nodes2[1]));
		REQUIRE(checkReprSet(nodes1, reprs1, 2));
		REQUIRE(checkReprSet(nodes2, reprs2, 2));

		// Tracing the distribution of MCC sizes by deleting edges one by one.
		DecrementalMCCAlgorithm::disconnect(edges1[0], reprs1, reprs2);
		REQUIRE(! SpanningForest::hasPath(nodes1[0], nodes1[1]));
		REQUIRE(! SpanningForest::hasPath(nodes2[0], nodes2[1]));
		REQUIRE(checkReprSet(nodes1, reprs1, 2));
		REQUIRE(checkReprSet(nodes2, reprs2, 2));
	}

	SECTION("") {
		// Configuration of network #1. Network #2 has no edges.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);

		// Initial distribution of MCC sizes.
		ClusterSizeDist cdist1;
		ClusterRepresentativeSet reprs1;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, cdist1, reprs1);
		REQUIRE(cdist1 == (ClusterSizeDist{{1, 2}}));
		REQUIRE(! SpanningForest::hasPath(nodes1[0], nodes1[1]));
		REQUIRE(! SpanningForest::hasPath(nodes2[0], nodes2[1]));
		REQUIRE(checkReprSet(nodes1, reprs1, 2));

		// Tracing the distribution of MCC sizes by deleting edges one by one.
		DecrementalMCCAlgorithm::disconnect(edges1[0], cdist1, reprs1);
		REQUIRE(cdist1 == (ClusterSizeDist{{1, 2}}));
		REQUIRE(! SpanningForest::hasPath(nodes1[0], nodes1[1]));
		REQUIRE(! SpanningForest::hasPath(nodes2[0], nodes2[1]));
		REQUIRE(checkReprSet(nodes1, reprs1, 2));
	}

	SECTION("") {
		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);

		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[1], edges2[0]);

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(SpanningForest::hasPath(nodes1[0], nodes1[1]));
		REQUIRE(SpanningForest::hasPath(nodes2[0], nodes2[1]));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));

		// Tracing the distribution of MCC sizes by deleting edges one by one.
		DecrementalMCCAlgorithm::disconnect(edges2[0], reprs1, reprs2);
		REQUIRE(! SpanningForest::hasPath(nodes1[0], nodes1[1]));
		REQUIRE(! SpanningForest::hasPath(nodes2[0], nodes2[1]));
		REQUIRE(checkReprSet(nodes1, reprs1, 2));
		REQUIRE(checkReprSet(nodes2, reprs2, 2));
	}

	SECTION("") {
		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);

		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[1], edges2[0]);

		// Initial distribution of MCC sizes.
		ClusterSizeDist cdist1;
		ClusterRepresentativeSet reprs1;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, cdist1, reprs1);
		REQUIRE(cdist1 == (ClusterSizeDist{{2, 1}}));
		REQUIRE(SpanningForest::hasPath(nodes1[0], nodes1[1]));
		REQUIRE(SpanningForest::hasPath(nodes2[0], nodes2[1]));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));

		// Tracing the distribution of MCC sizes by deleting edges one by one.
		DecrementalMCCAlgorithm::disconnect(edges2[0], cdist1, reprs1);
		REQUIRE(cdist1 == (ClusterSizeDist{{1, 2}}));
		REQUIRE(! SpanningForest::hasPath(nodes1[0], nodes1[1]));
		REQUIRE(! SpanningForest::hasPath(nodes2[0], nodes2[1]));
		REQUIRE(checkReprSet(nodes1, reprs1, 2));
	}
}

TEST_CASE("TestTwoNodesInterdependentNetworks-AllowMultipleInterdepedency") {
	{
		NodeVector nodes1(2); // nodes of network #1
		NodeVector nodes2(2); // nodes of network #2
		EdgeVector edges1(1); // edges of network #1
		EdgeVector edges2(1); // edges of network #2

		// Network #1 has no edges;
		// Network #2 has no edges;
		
		// Set interconnections between the two networks.
		// Each node has one or more interconnections.
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[0]);
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[1]);
		DecrementalMCCAlgorithm::interconnect(nodes1[1], nodes2[1]);

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(checkReprSet(nodes1, reprs1, 0));
		REQUIRE(checkReprSet(nodes2, reprs2, 0));
	}

	{
		NodeVector nodes1(2); // nodes of network #1
		NodeVector nodes2(2); // nodes of network #2
		EdgeVector edges1(1); // edges of network #1
		EdgeVector edges2(1); // edges of network #2

		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);

		// Network #2 has no edges;
		
		// Set interconnections between the two networks.
		// Each node has one or more interconnections.
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[0]);
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[1]);
		DecrementalMCCAlgorithm::interconnect(nodes1[1], nodes2[1]);

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(checkReprSet(nodes1, reprs1, 0));
		REQUIRE(checkReprSet(nodes2, reprs2, 0));
	}

	{
		NodeVector nodes1(2); // nodes of network #1
		NodeVector nodes2(2); // nodes of network #2
		EdgeVector edges1(1); // edges of network #1
		EdgeVector edges2(1); // edges of network #2

		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);

		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[1], edges2[0]);
		
		// Set interconnections between the two networks.
		// Each node has one or more interconnections.
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[0]);
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[1]);
		DecrementalMCCAlgorithm::interconnect(nodes1[1], nodes2[1]);

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));
	}
}	

TEST_CASE("TestThreeNodesInterdependentNetworks-SingleInterdepedency") {

	NodeVector nodes1(3); // nodes of network #1
	NodeVector nodes2(3); // nodes of network #2
	EdgeVector edges1(3); // edges of network #1
	EdgeVector edges2(3); // edges of network #2

	// Set interconnections between the two networks.
	// Each node has one interconnection.
	DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[0]);
	DecrementalMCCAlgorithm::interconnect(nodes1[1], nodes2[1]);
	DecrementalMCCAlgorithm::interconnect(nodes1[2], nodes2[2]);

	SECTION("") {
		// Each network contains no edges at all.

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(checkReprSet(nodes1, reprs1, 3));
		REQUIRE(checkReprSet(nodes2, reprs2, 3));
	}

	SECTION("") {
		// Each network contains no edges at all.

		// Initial distribution of MCC sizes.
		ClusterSizeDist cdist1;
		ClusterRepresentativeSet reprs1;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, cdist1, reprs1);
		REQUIRE(cdist1 == (ClusterSizeDist{{1, 3}}));
		REQUIRE(checkReprSet(nodes1, reprs1, 3));
	}

	SECTION("") {
		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);

		// Configuration of network #2.
		SpanningForest::connect(nodes2[1], nodes2[2], edges2[0]);

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 3));
		REQUIRE(checkReprSet(nodes2, reprs2, 3));
	}

	SECTION("") {
		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);

		// Configuration of network #2.
		SpanningForest::connect(nodes2[1], nodes2[2], edges2[0]);

		// Initial distribution of MCC sizes.
		ClusterSizeDist cdist1;
		ClusterRepresentativeSet reprs1;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, cdist1, reprs1);
		REQUIRE(cdist1 == (ClusterSizeDist{{1, 3}}));
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 3));
	}

	SECTION("") {
		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);

		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[1], edges2[0]);

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(SpanningForest::hasPath(nodes1[0], nodes1[1]));
		REQUIRE(! SpanningForest::hasPath(nodes1[1], nodes1[2]));
		REQUIRE(! SpanningForest::hasPath(nodes1[2], nodes1[0]));
		REQUIRE(SpanningForest::hasPath(nodes2[0], nodes2[1]));
		REQUIRE(! SpanningForest::hasPath(nodes2[1], nodes2[2]));
		REQUIRE(! SpanningForest::hasPath(nodes2[2], nodes2[0]));
		REQUIRE(checkReprSet(nodes1, reprs1, 2));
		REQUIRE(checkReprSet(nodes2, reprs2, 2));

		// Tracing the distribution of MCC sizes by deleting edges one by one.
		DecrementalMCCAlgorithm::disconnect(edges1[0], reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 3));
		REQUIRE(checkReprSet(nodes2, reprs2, 3));
	}

	SECTION("") {
		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);

		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[1], edges2[0]);

		// Initial distribution of MCC sizes.
		ClusterSizeDist cdist1;
		ClusterRepresentativeSet reprs1;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, cdist1, reprs1);
		REQUIRE(cdist1 == (ClusterSizeDist{{1, 1}, {2, 1}}));
		REQUIRE(SpanningForest::hasPath(nodes1[0], nodes1[1]));
		REQUIRE(! SpanningForest::hasPath(nodes1[1], nodes1[2]));
		REQUIRE(! SpanningForest::hasPath(nodes1[2], nodes1[0]));
		REQUIRE(SpanningForest::hasPath(nodes2[0], nodes2[1]));
		REQUIRE(! SpanningForest::hasPath(nodes2[1], nodes2[2]));
		REQUIRE(! SpanningForest::hasPath(nodes2[2], nodes2[0]));
		REQUIRE(checkReprSet(nodes1, reprs1, 2));

		// Tracing the distribution of MCC sizes by deleting edges one by one.
		DecrementalMCCAlgorithm::disconnect(edges1[0], cdist1, reprs1);
		REQUIRE(cdist1 == (ClusterSizeDist{{1, 3}}));
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 3));
	}

	SECTION("") {
		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);
		SpanningForest::connect(nodes1[1], nodes1[2], edges1[1]);

		// Configuration of network #2.
		SpanningForest::connect(nodes2[1], nodes2[2], edges2[0]);

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(! SpanningForest::hasPath(nodes1[0], nodes1[1]));
		REQUIRE(SpanningForest::hasPath(nodes1[1], nodes1[2]));
		REQUIRE(! SpanningForest::hasPath(nodes1[2], nodes1[0]));
		REQUIRE(! SpanningForest::hasPath(nodes2[0], nodes2[1]));
		REQUIRE(SpanningForest::hasPath(nodes2[1], nodes2[2]));
		REQUIRE(! SpanningForest::hasPath(nodes2[2], nodes2[0]));
		REQUIRE(checkReprSet(nodes1, reprs1, 2));
		REQUIRE(checkReprSet(nodes2, reprs2, 2));
	}

	SECTION("") {
		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);
		SpanningForest::connect(nodes1[1], nodes1[2], edges1[1]);

		// Configuration of network #2.
		SpanningForest::connect(nodes2[1], nodes2[2], edges2[0]);

		// Initial distribution of MCC sizes.
		ClusterSizeDist cdist1;
		ClusterRepresentativeSet reprs1;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, cdist1, reprs1);
		REQUIRE(cdist1 == (ClusterSizeDist{{1, 1}, {2, 1}}));
		REQUIRE(! SpanningForest::hasPath(nodes1[0], nodes1[1]));
		REQUIRE(SpanningForest::hasPath(nodes1[1], nodes1[2]));
		REQUIRE(! SpanningForest::hasPath(nodes1[2], nodes1[0]));
		REQUIRE(! SpanningForest::hasPath(nodes2[0], nodes2[1]));
		REQUIRE(SpanningForest::hasPath(nodes2[1], nodes2[2]));
		REQUIRE(! SpanningForest::hasPath(nodes2[2], nodes2[0]));
		REQUIRE(checkReprSet(nodes1, reprs1, 2));
	}

	SECTION("") {
		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);
		SpanningForest::connect(nodes1[1], nodes1[2], edges1[1]);
		SpanningForest::connect(nodes1[2], nodes1[0], edges1[2]);
		
		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[1], edges2[0]);
		SpanningForest::connect(nodes2[1], nodes2[2], edges2[1]);
		SpanningForest::connect(nodes2[2], nodes2[0], edges2[2]);
		
		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1, 2));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1, 2));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));

		// Tracing the distribution of MCC sizes by deleting edges one by one.
		DecrementalMCCAlgorithm::disconnect(edges1[0], reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1, 2));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1, 2));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));

		DecrementalMCCAlgorithm::disconnect(edges2[1], reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1, 2));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1, 2));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));

		DecrementalMCCAlgorithm::disconnect(edges1[1], reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 2));
		REQUIRE(cluster(nodes1[1]) == nodeSet(nodes1, 1));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 2));
		REQUIRE(cluster(nodes2[1]) == nodeSet(nodes2, 1));
		REQUIRE(checkReprSet(nodes1, reprs1, 2));
		REQUIRE(checkReprSet(nodes2, reprs2, 2));

		DecrementalMCCAlgorithm::disconnect(edges1[2], reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 3));
		REQUIRE(checkReprSet(nodes2, reprs2, 3));
	}

	SECTION("") {
		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);
		SpanningForest::connect(nodes1[1], nodes1[2], edges1[1]);
		SpanningForest::connect(nodes1[2], nodes1[0], edges1[2]);
		
		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[1], edges2[0]);
		SpanningForest::connect(nodes2[1], nodes2[2], edges2[1]);
		SpanningForest::connect(nodes2[2], nodes2[0], edges2[2]);
		
		// Initial distribution of MCC sizes.
		ClusterSizeDist cdist1;
		ClusterRepresentativeSet reprs1;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, cdist1, reprs1);
		REQUIRE(cdist1 == (ClusterSizeDist{{3, 1}}));
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1, 2));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1, 2));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));

		// Tracing the distribution of MCC sizes by deleting edges one by one.
		DecrementalMCCAlgorithm::disconnect(edges1[0], cdist1, reprs1);
		REQUIRE(cdist1 == (ClusterSizeDist{{3, 1}}));
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1, 2));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1, 2));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));

		DecrementalMCCAlgorithm::disconnect(edges2[1], cdist1, reprs1);
		REQUIRE(cdist1 == (ClusterSizeDist{{3, 1}}));
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1, 2));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1, 2));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));

		DecrementalMCCAlgorithm::disconnect(edges1[1], cdist1, reprs1);
		REQUIRE(cdist1 == (ClusterSizeDist{{1, 1}, {2, 1}}));
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 2));
		REQUIRE(cluster(nodes1[1]) == nodeSet(nodes1, 1));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 2));
		REQUIRE(cluster(nodes2[1]) == nodeSet(nodes2, 1));
		REQUIRE(checkReprSet(nodes1, reprs1, 2));

		DecrementalMCCAlgorithm::disconnect(edges1[2], cdist1, reprs1);
		REQUIRE(cdist1 == (ClusterSizeDist{{1, 3}}));
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 3));
	}
}

TEST_CASE("TestThreeNodesInterdependentNetworks-AllowMultipleInterdepedency") {

	{
		NodeVector nodes1(3); // nodes of network #1
		NodeVector nodes2(3); // nodes of network #2
		EdgeVector edges1(3); // edges of network #1
		EdgeVector edges2(3); // edges of network #2

		// Set interconnections between the two networks.
		// Each node has one or more interconnections.
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[0]);
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[1]);
		DecrementalMCCAlgorithm::interconnect(nodes1[1], nodes2[1]);
		DecrementalMCCAlgorithm::interconnect(nodes1[2], nodes2[2]);
	
		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);

		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[1], edges2[0]);

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1));
		REQUIRE(cluster(nodes1[2]) == nodeSet(nodes1, 2));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1));
		REQUIRE(cluster(nodes2[2]) == nodeSet(nodes2, 2));
		REQUIRE(checkReprSet(nodes1, reprs1, 2));
		REQUIRE(checkReprSet(nodes2, reprs2, 2));
	}

	{
		NodeVector nodes1(3); // nodes of network #1
		NodeVector nodes2(3); // nodes of network #2
		EdgeVector edges1(3); // edges of network #1
		EdgeVector edges2(3); // edges of network #2

		// Set interconnections between the two networks.
		// Each node has one or more interconnections.
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[0]);
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[2]);
		DecrementalMCCAlgorithm::interconnect(nodes1[1], nodes2[1]);
		DecrementalMCCAlgorithm::interconnect(nodes1[2], nodes2[2]);
	
		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);

		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[1], edges2[0]);

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));
	}
}

TEST_CASE("TestFourNodesInterdependentNetworks-SingleInterdepedency") {

	NodeVector nodes1(4); // nodes of network #1
	NodeVector nodes2(4); // nodes of network #2
	EdgeVector edges1(6); // edges of network #1
	EdgeVector edges2(6); // edges of network #2

	// Set interconnections between the two networks.
	// Each node has one interconnection.
	DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[0]);
	DecrementalMCCAlgorithm::interconnect(nodes1[1], nodes2[1]);
	DecrementalMCCAlgorithm::interconnect(nodes1[2], nodes2[2]);
	DecrementalMCCAlgorithm::interconnect(nodes1[3], nodes2[3]);

	SECTION("") {
		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);
		SpanningForest::connect(nodes1[1], nodes1[2], edges1[1]);
		SpanningForest::connect(nodes1[2], nodes1[0], edges1[2]);
		
		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[3], edges2[0]);
		SpanningForest::connect(nodes2[1], nodes2[3], edges2[1]);
		SpanningForest::connect(nodes2[2], nodes2[3], edges2[2]);
		
		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);

		// Tracing the distribution of MCC sizes by deleting edges one by one.
		DecrementalMCCAlgorithm::disconnect(edges1[0], reprs1, reprs2);
		REQUIRE(checkReprSet(nodes1, reprs1, 4));
		REQUIRE(checkReprSet(nodes2, reprs2, 4));

		DecrementalMCCAlgorithm::disconnect(edges2[0], reprs1, reprs2);
		REQUIRE(checkReprSet(nodes1, reprs1, 4));
		REQUIRE(checkReprSet(nodes2, reprs2, 4));
	}

	SECTION("") {
		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);
		SpanningForest::connect(nodes1[1], nodes1[2], edges1[1]);
		SpanningForest::connect(nodes1[2], nodes1[0], edges1[2]);
		
		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[3], edges2[0]);
		SpanningForest::connect(nodes2[1], nodes2[3], edges2[1]);
		SpanningForest::connect(nodes2[2], nodes2[3], edges2[2]);
		
		// Initial distribution of MCC sizes.
		ClusterSizeDist cdist1;
		ClusterRepresentativeSet reprs1;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, cdist1, reprs1);
		REQUIRE(cdist1 == (ClusterSizeDist{{1, 4}}));

		// Tracing the distribution of MCC sizes by deleting edges one by one.
		DecrementalMCCAlgorithm::disconnect(edges1[0], cdist1, reprs1);
		REQUIRE(cdist1 == (ClusterSizeDist{{1, 4}}));
		REQUIRE(checkReprSet(nodes1, reprs1, 4));

		DecrementalMCCAlgorithm::disconnect(edges2[0], cdist1, reprs1);
		REQUIRE(cdist1 == (ClusterSizeDist{{1, 4}}));
		REQUIRE(checkReprSet(nodes1, reprs1, 4));
	}

	SECTION("") {
		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);
		SpanningForest::connect(nodes1[1], nodes1[2], edges1[1]);
		SpanningForest::connect(nodes1[2], nodes1[0], edges1[2]);
		SpanningForest::connect(nodes1[1], nodes1[3], edges1[3]);
		
		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[3], edges2[0]);
		SpanningForest::connect(nodes2[1], nodes2[3], edges2[1]);
		SpanningForest::connect(nodes2[2], nodes2[3], edges2[2]);
		
		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1, 2, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));

		// Tracing the distribution of MCC sizes by deleting edges one by one.
		DecrementalMCCAlgorithm::disconnect(edges1[1], reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1, 2, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));

		DecrementalMCCAlgorithm::disconnect(edges1[3], reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 4));
		REQUIRE(checkReprSet(nodes2, reprs2, 4));
	}

	SECTION("") {
		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);
		SpanningForest::connect(nodes1[1], nodes1[2], edges1[1]);
		SpanningForest::connect(nodes1[2], nodes1[0], edges1[2]);
		SpanningForest::connect(nodes1[1], nodes1[3], edges1[3]);
		
		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[3], edges2[0]);
		SpanningForest::connect(nodes2[1], nodes2[3], edges2[1]);
		SpanningForest::connect(nodes2[2], nodes2[3], edges2[2]);
		
		// Initial distribution of MCC sizes.
		ClusterSizeDist cdist1;
		ClusterRepresentativeSet reprs1;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, cdist1, reprs1);
		REQUIRE(cdist1 == (ClusterSizeDist{{4, 1}}));
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1, 2, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));

		// Tracing the distribution of MCC sizes by deleting edges one by one.
		DecrementalMCCAlgorithm::disconnect(edges1[1], cdist1, reprs1);
		REQUIRE(cdist1 == (ClusterSizeDist{{4, 1}}));
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1, 2, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));

		DecrementalMCCAlgorithm::disconnect(edges1[3], cdist1, reprs1);
		REQUIRE(cdist1 == (ClusterSizeDist{{1, 4}}));
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 4));
	}

	SECTION("") {
		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);
		SpanningForest::connect(nodes1[1], nodes1[2], edges1[1]);
		SpanningForest::connect(nodes1[2], nodes1[0], edges1[2]);
		SpanningForest::connect(nodes1[1], nodes1[3], edges1[3]);
		
		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[3], edges2[0]);
		SpanningForest::connect(nodes2[1], nodes2[3], edges2[1]);
		SpanningForest::connect(nodes2[2], nodes2[3], edges2[2]);
		
		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1, 2, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));

		// Tracing the distribution of MCC sizes by deleting edges one by one.
		DecrementalMCCAlgorithm::disconnect(edges2[0], reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0));
		REQUIRE(cluster(nodes1[1]) == nodeSet(nodes1, 1, 2, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0));
		REQUIRE(cluster(nodes2[1]) == nodeSet(nodes2, 1, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 2));
		REQUIRE(checkReprSet(nodes2, reprs2, 2));

		DecrementalMCCAlgorithm::disconnect(edges1[1], reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0));
		REQUIRE(cluster(nodes1[1]) == nodeSet(nodes1, 1, 3));
		REQUIRE(cluster(nodes1[2]) == nodeSet(nodes1, 2));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0));
		REQUIRE(cluster(nodes2[1]) == nodeSet(nodes2, 1, 3));
		REQUIRE(cluster(nodes2[2]) == nodeSet(nodes2, 2));
		REQUIRE(checkReprSet(nodes1, reprs1, 3));
		REQUIRE(checkReprSet(nodes2, reprs2, 3));

		DecrementalMCCAlgorithm::disconnect(edges2[1], reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 4));
		REQUIRE(checkReprSet(nodes2, reprs2, 4));
	}

	SECTION("") {
		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);
		SpanningForest::connect(nodes1[1], nodes1[2], edges1[1]);
		SpanningForest::connect(nodes1[2], nodes1[0], edges1[2]);
		SpanningForest::connect(nodes1[1], nodes1[3], edges1[3]);
		
		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[3], edges2[0]);
		SpanningForest::connect(nodes2[1], nodes2[3], edges2[1]);
		SpanningForest::connect(nodes2[2], nodes2[3], edges2[2]);
		
		// Initial distribution of MCC sizes.
		ClusterSizeDist cdist1;
		ClusterRepresentativeSet reprs1;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, cdist1, reprs1);
		REQUIRE(cdist1 == (ClusterSizeDist{{4, 1}}));
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1, 2, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));

		// Tracing the distribution of MCC sizes by deleting edges one by one.
		DecrementalMCCAlgorithm::disconnect(edges2[0], cdist1, reprs1);
		REQUIRE(cdist1 == (ClusterSizeDist{{1, 1}, {3, 1}}));
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0));
		REQUIRE(cluster(nodes1[1]) == nodeSet(nodes1, 1, 2, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0));
		REQUIRE(cluster(nodes2[1]) == nodeSet(nodes2, 1, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 2));

		DecrementalMCCAlgorithm::disconnect(edges1[1], cdist1, reprs1);
		REQUIRE(cdist1 == (ClusterSizeDist{{1, 2}, {2, 1}}));
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0));
		REQUIRE(cluster(nodes1[1]) == nodeSet(nodes1, 1, 3));
		REQUIRE(cluster(nodes1[2]) == nodeSet(nodes1, 2));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0));
		REQUIRE(cluster(nodes2[1]) == nodeSet(nodes2, 1, 3));
		REQUIRE(cluster(nodes2[2]) == nodeSet(nodes2, 2));
		REQUIRE(checkReprSet(nodes1, reprs1, 3));

		DecrementalMCCAlgorithm::disconnect(edges2[1], cdist1, reprs1);
		REQUIRE(cdist1 == (ClusterSizeDist{{1, 4}}));
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 4));
	}

	SECTION("") {
		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);
		SpanningForest::connect(nodes1[0], nodes1[2], edges1[1]);
		SpanningForest::connect(nodes1[1], nodes1[2], edges1[2]);
		SpanningForest::connect(nodes1[1], nodes1[3], edges1[3]);
		SpanningForest::connect(nodes1[2], nodes1[3], edges1[4]);
		
		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[1], edges2[0]);
		SpanningForest::connect(nodes2[0], nodes2[2], edges2[1]);
		SpanningForest::connect(nodes2[1], nodes2[2], edges2[2]);
		SpanningForest::connect(nodes2[2], nodes2[3], edges2[3]);
		
		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1, 2, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));

		// Tracing the distribution of MCC sizes by deleting edges one by one.
		DecrementalMCCAlgorithm::disconnect(edges2[1], reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1, 2, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));

		DecrementalMCCAlgorithm::disconnect(edges2[2], reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1));
		REQUIRE(cluster(nodes1[2]) == nodeSet(nodes1, 2, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1));
		REQUIRE(cluster(nodes2[2]) == nodeSet(nodes2, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 2));
		REQUIRE(checkReprSet(nodes2, reprs2, 2));

		DecrementalMCCAlgorithm::disconnect(edges1[1], reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1));
		REQUIRE(cluster(nodes1[2]) == nodeSet(nodes1, 2, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1));
		REQUIRE(cluster(nodes2[2]) == nodeSet(nodes2, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 2));
		REQUIRE(checkReprSet(nodes2, reprs2, 2));

		DecrementalMCCAlgorithm::disconnect(edges1[0], reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0));
		REQUIRE(cluster(nodes1[1]) == nodeSet(nodes1, 1));
		REQUIRE(cluster(nodes1[2]) == nodeSet(nodes1, 2, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0));
		REQUIRE(cluster(nodes2[1]) == nodeSet(nodes2, 1));
		REQUIRE(cluster(nodes2[2]) == nodeSet(nodes2, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 3));
		REQUIRE(checkReprSet(nodes2, reprs2, 3));

		DecrementalMCCAlgorithm::disconnect(edges1[2], reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0));
		REQUIRE(cluster(nodes1[1]) == nodeSet(nodes1, 1));
		REQUIRE(cluster(nodes1[2]) == nodeSet(nodes1, 2, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0));
		REQUIRE(cluster(nodes2[1]) == nodeSet(nodes2, 1));
		REQUIRE(cluster(nodes2[2]) == nodeSet(nodes2, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 3));
		REQUIRE(checkReprSet(nodes2, reprs2, 3));

		DecrementalMCCAlgorithm::disconnect(edges1[3], reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0));
		REQUIRE(cluster(nodes1[1]) == nodeSet(nodes1, 1));
		REQUIRE(cluster(nodes1[2]) == nodeSet(nodes1, 2, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0));
		REQUIRE(cluster(nodes2[1]) == nodeSet(nodes2, 1));
		REQUIRE(cluster(nodes2[2]) == nodeSet(nodes2, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 3));
		REQUIRE(checkReprSet(nodes2, reprs2, 3));

		DecrementalMCCAlgorithm::disconnect(edges1[4], reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 4));
		REQUIRE(checkReprSet(nodes2, reprs2, 4));
	}

	SECTION("") {
		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);
		SpanningForest::connect(nodes1[0], nodes1[2], edges1[1]);
		SpanningForest::connect(nodes1[1], nodes1[2], edges1[2]);
		SpanningForest::connect(nodes1[1], nodes1[3], edges1[3]);
		SpanningForest::connect(nodes1[2], nodes1[3], edges1[4]);
		
		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[1], edges2[0]);
		SpanningForest::connect(nodes2[0], nodes2[2], edges2[1]);
		SpanningForest::connect(nodes2[1], nodes2[2], edges2[2]);
		SpanningForest::connect(nodes2[2], nodes2[3], edges2[3]);
		
		// Initial distribution of MCC sizes.
		ClusterSizeDist cdist1;
		ClusterRepresentativeSet reprs1;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, cdist1, reprs1);
		REQUIRE(cdist1 == (ClusterSizeDist{{4, 1}}));
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1, 2, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));

		// Tracing the distribution of MCC sizes by deleting edges one by one.
		DecrementalMCCAlgorithm::disconnect(edges2[1], cdist1, reprs1);
		REQUIRE(cdist1 == (ClusterSizeDist{{4, 1}}));
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1, 2, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));

		DecrementalMCCAlgorithm::disconnect(edges2[2], cdist1, reprs1);
		REQUIRE(cdist1 == (ClusterSizeDist{{2, 2}}));
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1));
		REQUIRE(cluster(nodes1[2]) == nodeSet(nodes1, 2, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1));
		REQUIRE(cluster(nodes2[2]) == nodeSet(nodes2, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 2));

		DecrementalMCCAlgorithm::disconnect(edges1[1], cdist1, reprs1);
		REQUIRE(cdist1 == (ClusterSizeDist{{2, 2}}));
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1));
		REQUIRE(cluster(nodes1[2]) == nodeSet(nodes1, 2, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1));
		REQUIRE(cluster(nodes2[2]) == nodeSet(nodes2, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 2));

		DecrementalMCCAlgorithm::disconnect(edges1[0], cdist1, reprs1);
		REQUIRE(cdist1 == (ClusterSizeDist{{1, 2}, {2, 1}}));
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0));
		REQUIRE(cluster(nodes1[1]) == nodeSet(nodes1, 1));
		REQUIRE(cluster(nodes1[2]) == nodeSet(nodes1, 2, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0));
		REQUIRE(cluster(nodes2[1]) == nodeSet(nodes2, 1));
		REQUIRE(cluster(nodes2[2]) == nodeSet(nodes2, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 3));

		DecrementalMCCAlgorithm::disconnect(edges1[2], cdist1, reprs1);
		REQUIRE(cdist1 == (ClusterSizeDist{{1, 2}, {2, 1}}));
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0));
		REQUIRE(cluster(nodes1[1]) == nodeSet(nodes1, 1));
		REQUIRE(cluster(nodes1[2]) == nodeSet(nodes1, 2, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0));
		REQUIRE(cluster(nodes2[1]) == nodeSet(nodes2, 1));
		REQUIRE(cluster(nodes2[2]) == nodeSet(nodes2, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 3));

		DecrementalMCCAlgorithm::disconnect(edges1[3], cdist1, reprs1);
		REQUIRE(cdist1 == (ClusterSizeDist{{1, 2}, {2, 1}}));
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0));
		REQUIRE(cluster(nodes1[1]) == nodeSet(nodes1, 1));
		REQUIRE(cluster(nodes1[2]) == nodeSet(nodes1, 2, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0));
		REQUIRE(cluster(nodes2[1]) == nodeSet(nodes2, 1));
		REQUIRE(cluster(nodes2[2]) == nodeSet(nodes2, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 3));

		DecrementalMCCAlgorithm::disconnect(edges1[4], cdist1, reprs1);
		REQUIRE(cdist1 == (ClusterSizeDist{{1, 4}}));
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 4));
	}

	SECTION("") {
		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);
		SpanningForest::connect(nodes1[0], nodes1[2], edges1[1]);
		SpanningForest::connect(nodes1[0], nodes1[3], edges1[2]);
		
		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[1], edges2[0]);
		SpanningForest::connect(nodes2[1], nodes2[2], edges2[1]);
		SpanningForest::connect(nodes2[1], nodes2[3], edges2[2]);
		SpanningForest::connect(nodes2[2], nodes2[3], edges2[3]);
		
		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1, 2, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));

		// Tracing the distribution of MCC sizes by deleting edges one by one.
		DecrementalMCCAlgorithm::disconnect(nodes1[0], reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 4));
		REQUIRE(checkReprSet(nodes2, reprs2, 4));
	}

	SECTION("") {
		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);
		SpanningForest::connect(nodes1[0], nodes1[2], edges1[1]);
		SpanningForest::connect(nodes1[0], nodes1[3], edges1[2]);
		
		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[1], edges2[0]);
		SpanningForest::connect(nodes2[1], nodes2[2], edges2[1]);
		SpanningForest::connect(nodes2[1], nodes2[3], edges2[2]);
		SpanningForest::connect(nodes2[2], nodes2[3], edges2[3]);
		
		// Initial distribution of MCC sizes.
		ClusterSizeDist cdist1;
		ClusterRepresentativeSet reprs1;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, cdist1, reprs1);
		REQUIRE(cdist1 == (ClusterSizeDist{{4, 1}}));
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1, 2, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));

		// Tracing the distribution of MCC sizes by deleting edges one by one.
		DecrementalMCCAlgorithm::disconnect(nodes1[0], cdist1, reprs1);
		REQUIRE(cdist1 == (ClusterSizeDist{{1, 4}}));
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 4));
	}
}

TEST_CASE("TestFourNodesInterdependentNetworks-AllowMultipleInterdepedency") {

	{
		NodeVector nodes1(4); // nodes of network #1
		NodeVector nodes2(4); // nodes of network #2
		EdgeVector edges1(6); // edges of network #1
		EdgeVector edges2(6); // edges of network #2

		// Configuration of network #1.
		SpanningForest::connect(nodes1[1], nodes1[2], edges1[0]);
		SpanningForest::connect(nodes1[2], nodes1[3], edges1[1]);
		
		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[1], edges2[0]);
		SpanningForest::connect(nodes2[1], nodes2[2], edges2[1]);
		SpanningForest::connect(nodes2[2], nodes2[3], edges2[2]);

		// Set interconnections between the two networks.
		// Each node has one or more interconnections.
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[0]);
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[1]);
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[2]);
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[3]);
		DecrementalMCCAlgorithm::interconnect(nodes1[1], nodes2[1]);
		DecrementalMCCAlgorithm::interconnect(nodes1[2], nodes2[2]);
		DecrementalMCCAlgorithm::interconnect(nodes1[3], nodes2[3]);

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 0));
		REQUIRE(checkReprSet(nodes2, reprs2, 0));
	}

	{
		NodeVector nodes1(4); // nodes of network #1
		NodeVector nodes2(4); // nodes of network #2
		EdgeVector edges1(6); // edges of network #1
		EdgeVector edges2(6); // edges of network #2

		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);
		SpanningForest::connect(nodes1[1], nodes1[2], edges1[1]);
		SpanningForest::connect(nodes1[2], nodes1[3], edges1[2]);
		
		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[1], edges2[0]);
		SpanningForest::connect(nodes2[1], nodes2[2], edges2[1]);
		SpanningForest::connect(nodes2[2], nodes2[3], edges2[2]);

		// Set interconnections between the two networks.
		// Each node has one or more interconnections.
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[0]);
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[1]);
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[2]);
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[3]);
		DecrementalMCCAlgorithm::interconnect(nodes1[1], nodes2[1]);
		DecrementalMCCAlgorithm::interconnect(nodes1[2], nodes2[2]);
		DecrementalMCCAlgorithm::interconnect(nodes1[3], nodes2[3]);

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1, 2, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));

		// Tracing the distribution of MCC sizes by deleting edges one by one.
		DecrementalMCCAlgorithm::disconnect(edges1[0], reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 0));
		REQUIRE(checkReprSet(nodes2, reprs2, 0));
	}

	{
		NodeVector nodes1(4); // nodes of network #1
		NodeVector nodes2(4); // nodes of network #2
		EdgeVector edges1(6); // edges of network #1
		EdgeVector edges2(6); // edges of network #2

		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);
		SpanningForest::connect(nodes1[1], nodes1[2], edges1[1]);
		SpanningForest::connect(nodes1[2], nodes1[3], edges1[2]);
		SpanningForest::connect(nodes1[3], nodes1[1], edges1[3]);
		
		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[1], edges2[0]);
		SpanningForest::connect(nodes2[1], nodes2[2], edges2[1]);
		SpanningForest::connect(nodes2[2], nodes2[3], edges2[2]);
		SpanningForest::connect(nodes2[3], nodes2[1], edges2[3]);

		// Set interconnections between the two networks.
		// Each node has one or more interconnections.
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[0]);
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[1]);
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[2]);
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[3]);
		DecrementalMCCAlgorithm::interconnect(nodes1[1], nodes2[1]);
		DecrementalMCCAlgorithm::interconnect(nodes1[2], nodes2[2]);
		DecrementalMCCAlgorithm::interconnect(nodes1[3], nodes2[3]);

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1, 2, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));

		// Tracing the distribution of MCC sizes by deleting edges one by one.
		DecrementalMCCAlgorithm::disconnect(edges1[1], reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1, 2, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));

		DecrementalMCCAlgorithm::disconnect(edges2[1], reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1, 2, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));

		DecrementalMCCAlgorithm::disconnect(edges1[0], reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 0));
		REQUIRE(checkReprSet(nodes2, reprs2, 0));
	}

	{
		NodeVector nodes1(3); // nodes of network #1
		NodeVector nodes2(3); // nodes of network #2
		EdgeVector edges1(3); // edges of network #1
		EdgeVector edges2(3); // edges of network #2

		// Configuration of network #1.
		SpanningForest::connect(nodes1[1], nodes1[2], edges1[0]);
		
		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[1], edges2[0]);

		// Set interconnections between the two networks.
		// Each node has one or more interconnections.
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[0]);
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[1]);
		DecrementalMCCAlgorithm::interconnect(nodes1[1], nodes2[2]);
		DecrementalMCCAlgorithm::interconnect(nodes1[2], nodes2[2]);

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0));
		REQUIRE(cluster(nodes1[1]) == nodeSet(nodes1, 1, 2));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1));
		REQUIRE(cluster(nodes2[2]) == nodeSet(nodes2, 2));
		REQUIRE(checkReprSet(nodes1, reprs1, 2));
		REQUIRE(checkReprSet(nodes2, reprs2, 2));

		// Tracing the distribution of MCC sizes by deleting edges one by one.
		DecrementalMCCAlgorithm::disconnect(edges1[0], reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0));
		REQUIRE(cluster(nodes1[1]) == nodeSet(nodes1, 1));
		REQUIRE(cluster(nodes1[2]) == nodeSet(nodes1, 2));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1));
		REQUIRE(cluster(nodes2[2]) == nodeSet(nodes2, 2));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));

		DecrementalMCCAlgorithm::disconnect(edges2[0], reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 0));
		REQUIRE(checkReprSet(nodes2, reprs2, 0));
	}

	{
		NodeVector nodes1(4); // nodes of network #1
		NodeVector nodes2(4); // nodes of network #2
		EdgeVector edges1(6); // edges of network #1
		EdgeVector edges2(6); // edges of network #2

		// Configuration of network #1.
		SpanningForest::connect(nodes1[1], nodes1[2], edges1[0]);
		SpanningForest::connect(nodes1[2], nodes1[3], edges1[1]);
		
		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[1], edges2[0]);
		SpanningForest::connect(nodes2[1], nodes2[2], edges2[1]);

		// Set interconnections between the two networks.
		// Each node has one or more interconnections.
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[0]);
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[1]);
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[2]);
		DecrementalMCCAlgorithm::interconnect(nodes1[1], nodes2[3]);
		DecrementalMCCAlgorithm::interconnect(nodes1[2], nodes2[3]);
		DecrementalMCCAlgorithm::interconnect(nodes1[3], nodes2[3]);

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0));
		REQUIRE(cluster(nodes1[1]) == nodeSet(nodes1, 1, 2, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1, 2));
		REQUIRE(cluster(nodes2[3]) == nodeSet(nodes2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 2));
		REQUIRE(checkReprSet(nodes2, reprs2, 2));

		// Tracing the distribution of MCC sizes by deleting edges one by one.
		DecrementalMCCAlgorithm::disconnect(edges1[0], reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0));
		REQUIRE(cluster(nodes1[1]) == nodeSet(nodes1, 1));
		REQUIRE(cluster(nodes1[2]) == nodeSet(nodes1, 2));
		REQUIRE(cluster(nodes1[3]) == nodeSet(nodes1, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1, 2));
		REQUIRE(cluster(nodes2[3]) == nodeSet(nodes2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));

		DecrementalMCCAlgorithm::disconnect(edges2[1], reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 0));
		REQUIRE(checkReprSet(nodes2, reprs2, 0));
	}

	{
		NodeVector nodes1(4); // nodes of network #1
		NodeVector nodes2(4); // nodes of network #2
		EdgeVector edges1(6); // edges of network #1
		EdgeVector edges2(6); // edges of network #2

		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);
		SpanningForest::connect(nodes1[1], nodes1[2], edges1[1]);
		SpanningForest::connect(nodes1[2], nodes1[0], edges1[2]);
		SpanningForest::connect(nodes1[1], nodes1[3], edges1[3]);
		SpanningForest::connect(nodes1[2], nodes1[3], edges1[4]);
		
		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[3], edges2[0]);
		SpanningForest::connect(nodes2[1], nodes2[3], edges2[1]);
		SpanningForest::connect(nodes2[2], nodes2[3], edges2[2]);

		// Set interconnections between the two networks.
		// Each node has one or more interconnections.
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[0]);
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[1]);
		DecrementalMCCAlgorithm::interconnect(nodes1[1], nodes2[1]);
		DecrementalMCCAlgorithm::interconnect(nodes1[2], nodes2[2]);
		DecrementalMCCAlgorithm::interconnect(nodes1[3], nodes2[3]);

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1, 2, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));

		// Tracing the distribution of MCC sizes by deleting edges one by one.
		DecrementalMCCAlgorithm::disconnect(edges2[0], reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0));
		REQUIRE(cluster(nodes1[1]) == nodeSet(nodes1, 1));
		REQUIRE(cluster(nodes1[2]) == nodeSet(nodes1, 2, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0));
		REQUIRE(cluster(nodes2[1]) == nodeSet(nodes2, 1));
		REQUIRE(cluster(nodes2[2]) == nodeSet(nodes2, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));
	}

	{
		NodeVector nodes1(4); // nodes of network #1
		NodeVector nodes2(4); // nodes of network #2
		EdgeVector edges1(6); // edges of network #1
		EdgeVector edges2(6); // edges of network #2

		// Set interconnections between the two networks.
		// Each node has one or more interconnections.
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[0]);
		DecrementalMCCAlgorithm::interconnect(nodes1[1], nodes2[1]);
		DecrementalMCCAlgorithm::interconnect(nodes1[1], nodes2[3]);
		DecrementalMCCAlgorithm::interconnect(nodes1[2], nodes2[2]);
		DecrementalMCCAlgorithm::interconnect(nodes1[3], nodes2[3]);

		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);
		SpanningForest::connect(nodes1[1], nodes1[2], edges1[1]);
		SpanningForest::connect(nodes1[2], nodes1[0], edges1[2]);
		SpanningForest::connect(nodes1[1], nodes1[3], edges1[3]);
		
		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[3], edges2[0]);
		SpanningForest::connect(nodes2[1], nodes2[3], edges2[1]);
		SpanningForest::connect(nodes2[2], nodes2[3], edges2[2]);

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1, 2, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));

		// Tracing the distribution of MCC sizes by deleting edges one by one.
		DecrementalMCCAlgorithm::disconnect(edges2[0], reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0));
		REQUIRE(cluster(nodes1[1]) == nodeSet(nodes1, 1, 2, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0));
		REQUIRE(cluster(nodes2[1]) == nodeSet(nodes2, 1, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 2));
		REQUIRE(checkReprSet(nodes2, reprs2, 2));

		DecrementalMCCAlgorithm::disconnect(edges1[1], reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0));
		REQUIRE(cluster(nodes1[1]) == nodeSet(nodes1, 1, 3));
		REQUIRE(cluster(nodes1[2]) == nodeSet(nodes1, 2));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0));
		REQUIRE(cluster(nodes2[1]) == nodeSet(nodes2, 1, 3));
		REQUIRE(cluster(nodes2[2]) == nodeSet(nodes2, 2));
		REQUIRE(checkReprSet(nodes1, reprs1, 3));
		REQUIRE(checkReprSet(nodes2, reprs2, 3));
	}

	{
		NodeVector nodes1(4); // nodes of network #1
		NodeVector nodes2(4); // nodes of network #2
		EdgeVector edges1(6); // edges of network #1
		EdgeVector edges2(6); // edges of network #2

		// Set interconnections between the two networks.
		// Each node has one or more interconnections.
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[0]);
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[1]);
		DecrementalMCCAlgorithm::interconnect(nodes1[1], nodes2[0]);
		DecrementalMCCAlgorithm::interconnect(nodes1[1], nodes2[1]);
		DecrementalMCCAlgorithm::interconnect(nodes1[2], nodes2[2]);
		DecrementalMCCAlgorithm::interconnect(nodes1[2], nodes2[3]);
		DecrementalMCCAlgorithm::interconnect(nodes1[3], nodes2[2]);
		DecrementalMCCAlgorithm::interconnect(nodes1[3], nodes2[3]);

		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);
		SpanningForest::connect(nodes1[1], nodes1[3], edges1[1]);
		SpanningForest::connect(nodes1[1], nodes1[2], edges1[2]);
		SpanningForest::connect(nodes1[2], nodes1[3], edges1[3]);
		
		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[1], edges2[0]);
		SpanningForest::connect(nodes2[2], nodes2[3], edges2[1]);

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1));
		REQUIRE(cluster(nodes1[2]) == nodeSet(nodes1, 2, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1));
		REQUIRE(cluster(nodes2[2]) == nodeSet(nodes2, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 2));
		REQUIRE(checkReprSet(nodes2, reprs2, 2));
	}

	{
		NodeVector nodes1(4); // nodes of network #1
		NodeVector nodes2(4); // nodes of network #2
		EdgeVector edges1(6); // edges of network #1
		EdgeVector edges2(6); // edges of network #2

		// Set interconnections between the two networks.
		// Each node has one or more interconnections.
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[0]);
		DecrementalMCCAlgorithm::interconnect(nodes1[1], nodes2[1]);
		DecrementalMCCAlgorithm::interconnect(nodes1[2], nodes2[2]);
		DecrementalMCCAlgorithm::interconnect(nodes1[3], nodes2[3]);
		DecrementalMCCAlgorithm::interconnect(nodes1[1], nodes2[2]);

		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);
		SpanningForest::connect(nodes1[2], nodes1[3], edges1[1]);
		
		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[1], edges2[0]);
		SpanningForest::connect(nodes2[2], nodes2[3], edges2[1]);

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 2));
		REQUIRE(checkReprSet(nodes2, reprs2, 2));

		// Tracing the distribution of MCC sizes by deleting edges one by one.
		DecrementalMCCAlgorithm::disconnect(edges1[0], reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 2));
		REQUIRE(checkReprSet(nodes2, reprs2, 2));

		DecrementalMCCAlgorithm::disconnect(edges1[1], reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 2));
		REQUIRE(checkReprSet(nodes2, reprs2, 2));

		DecrementalMCCAlgorithm::disconnect(edges2[0], reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 2));
		REQUIRE(checkReprSet(nodes2, reprs2, 2));

		DecrementalMCCAlgorithm::disconnect(edges2[1], reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 2));
		REQUIRE(checkReprSet(nodes2, reprs2, 2));
	}

	{
		NodeVector nodes1(4); // nodes of network #1
		NodeVector nodes2(4); // nodes of network #2
		EdgeVector edges1(6); // edges of network #1
		EdgeVector edges2(6); // edges of network #2

		// Set interconnections between the two networks.
		// Each node has one or more interconnections.
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[0]);
		DecrementalMCCAlgorithm::interconnect(nodes1[1], nodes2[1]);
		DecrementalMCCAlgorithm::interconnect(nodes1[2], nodes2[2]);
		DecrementalMCCAlgorithm::interconnect(nodes1[3], nodes2[3]);
		DecrementalMCCAlgorithm::interconnect(nodes1[1], nodes2[2]);
		DecrementalMCCAlgorithm::interconnect(nodes1[2], nodes2[3]);

		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);
		SpanningForest::connect(nodes1[2], nodes1[3], edges1[1]);
		
		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[1], edges2[0]);
		SpanningForest::connect(nodes2[2], nodes2[3], edges2[1]);

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));

		DecrementalMCCAlgorithm::disconnect(edges1[0], reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));

		DecrementalMCCAlgorithm::disconnect(edges1[1], reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));

		DecrementalMCCAlgorithm::disconnect(edges2[0], reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));

		DecrementalMCCAlgorithm::disconnect(edges2[1], reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));
	}
}

TEST_CASE("TestFiveNodesInterdependentNetworks-AllowMultipleInterdepedency") {

	{
		NodeVector nodes1(5); // nodes of network #1
		NodeVector nodes2(5); // nodes of network #2
		EdgeVector edges1(3); // edges of network #1
		EdgeVector edges2(3); // edges of network #2

		// Set interconnections between the two networks.
		// Each node has one or more interconnections.
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[0]);
		DecrementalMCCAlgorithm::interconnect(nodes1[1], nodes2[1]);
		DecrementalMCCAlgorithm::interconnect(nodes1[2], nodes2[2]);
		DecrementalMCCAlgorithm::interconnect(nodes1[3], nodes2[3]);
		DecrementalMCCAlgorithm::interconnect(nodes1[4], nodes2[4]);
		DecrementalMCCAlgorithm::interconnect(nodes1[1], nodes2[2]);
		DecrementalMCCAlgorithm::interconnect(nodes1[2], nodes2[3]);
		DecrementalMCCAlgorithm::interconnect(nodes1[3], nodes2[4]);

		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);
		SpanningForest::connect(nodes1[2], nodes1[3], edges1[1]);
		SpanningForest::connect(nodes1[3], nodes1[4], edges1[2]);
		
		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[1], edges2[0]);
		SpanningForest::connect(nodes2[2], nodes2[3], edges2[1]);
		SpanningForest::connect(nodes2[3], nodes2[4], edges2[2]);

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));

		DecrementalMCCAlgorithm::disconnect(edges1[0], reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));

		DecrementalMCCAlgorithm::disconnect(edges1[1], reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));

		DecrementalMCCAlgorithm::disconnect(edges2[0], reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));

		DecrementalMCCAlgorithm::disconnect(edges2[1], reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));
	}
}

TEST_CASE("TestInterdependentNetworks-AllowNoInterdepedencyAndDifferentNetworkSizes") {

	{
		// Trivial single node networks
		NodeVector nodes1(1);
		NodeVector nodes2(1);
		EdgeVector edges1;
		EdgeVector edges2;

		// No interconnections

		// Each network is totally disconnected.

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(checkReprSet(nodes1, reprs1, 0));
		REQUIRE(checkReprSet(nodes2, reprs2, 0));
	}

	{
		// Trivial single node networks
		NodeVector nodes1(1);
		NodeVector nodes2(2);
		EdgeVector edges1;
		EdgeVector edges2;

		// No interconnections

		// Each network is totally disconnected.

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(checkReprSet(nodes1, reprs1, 0));
		REQUIRE(checkReprSet(nodes2, reprs2, 0));
	}

	{
		NodeVector nodes1(2); // nodes of network #1
		NodeVector nodes2(2); // nodes of network #2
		EdgeVector edges1(1); // edges of network #1
		EdgeVector edges2(1); // edges of network #2

		// No interconnections

		// Each network is totally disconnected.

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(checkReprSet(nodes1, reprs1, 0));
		REQUIRE(checkReprSet(nodes2, reprs2, 0));
	}

	{
		NodeVector nodes1(2); // nodes of network #1
		NodeVector nodes2(2); // nodes of network #2
		EdgeVector edges1(1); // edges of network #1
		EdgeVector edges2(1); // edges of network #2

		// Set interconnections between the two networks.
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[0]);

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));
	}

	{
		NodeVector nodes1(1); // nodes of network #1
		NodeVector nodes2(2); // nodes of network #2
		EdgeVector edges1(0); // edges of network #1
		EdgeVector edges2(1); // edges of network #2

		// Set interconnections between the two networks.
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[0]);

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));
	}

	{
		NodeVector nodes1(2); // nodes of network #1
		NodeVector nodes2(1); // nodes of network #2
		EdgeVector edges1(0); // edges of network #1
		EdgeVector edges2(1); // edges of network #2

		// Set interconnections between the two networks.
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[0]);

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));
	}

	{
		NodeVector nodes1(1); // nodes of network #1
		NodeVector nodes2(2); // nodes of network #2
		EdgeVector edges1(0); // edges of network #1
		EdgeVector edges2(1); // edges of network #2

		// Set interconnections between the two networks.
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[0]);

		// Network #1 has no edges.

		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[1], edges2[0]);

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));
	}

	{
		NodeVector nodes1(4); // nodes of network #1
		NodeVector nodes2(2); // nodes of network #2
		EdgeVector edges1(6); // edges of network #1
		EdgeVector edges2(1); // edges of network #2

		// Set interconnections between the two networks.
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[0]);
		DecrementalMCCAlgorithm::interconnect(nodes1[1], nodes2[1]);

		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);
		SpanningForest::connect(nodes1[1], nodes1[2], edges1[1]);
		SpanningForest::connect(nodes1[2], nodes1[3], edges1[2]);

		// Network #2 has no edges.

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 2));
		REQUIRE(checkReprSet(nodes2, reprs2, 2));
	}

	{
		NodeVector nodes1(4); // nodes of network #1
		NodeVector nodes2(2); // nodes of network #2
		EdgeVector edges1(6); // edges of network #1
		EdgeVector edges2(1); // edges of network #2

		// Set interconnections between the two networks.
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[0]);
		DecrementalMCCAlgorithm::interconnect(nodes1[1], nodes2[1]);

		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);
		SpanningForest::connect(nodes1[1], nodes1[2], edges1[1]);
		SpanningForest::connect(nodes1[2], nodes1[3], edges1[2]);

		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[1], edges1[0]);

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1));
		REQUIRE(cluster(nodes1[2]) == nodeSet(nodes1, 2));
		REQUIRE(cluster(nodes1[3]) == nodeSet(nodes1, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));
	}

	{
		NodeVector nodes1(4); // nodes of network #1
		NodeVector nodes2(2); // nodes of network #2
		EdgeVector edges1(6); // edges of network #1
		EdgeVector edges2(1); // edges of network #2

		// Set interconnections between the two networks.
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[0]);
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[1]);
		DecrementalMCCAlgorithm::interconnect(nodes1[1], nodes2[1]);

		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);
		SpanningForest::connect(nodes1[1], nodes1[2], edges1[1]);
		SpanningForest::connect(nodes1[2], nodes1[3], edges1[2]);

		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[1], edges1[0]);

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1));
		REQUIRE(cluster(nodes1[2]) == nodeSet(nodes1, 2));
		REQUIRE(cluster(nodes1[3]) == nodeSet(nodes1, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));
	}

	{
		NodeVector nodes1(4); // nodes of network #1
		NodeVector nodes2(2); // nodes of network #2
		EdgeVector edges1(6); // edges of network #1
		EdgeVector edges2(1); // edges of network #2

		// Set interconnections between the two networks.
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[0]);
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[1]);
		DecrementalMCCAlgorithm::interconnect(nodes1[1], nodes2[1]);

		// Configuration of network #1.
		SpanningForest::connect(nodes1[1], nodes1[2], edges1[1]);
		SpanningForest::connect(nodes1[2], nodes1[3], edges1[2]);

		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[1], edges1[0]);

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 0));
		REQUIRE(checkReprSet(nodes2, reprs2, 0));
	}

	{
		NodeVector nodes1(3); // nodes of network #1
		NodeVector nodes2(3); // nodes of network #2
		EdgeVector edges1(3); // edges of network #1
		EdgeVector edges2(3); // edges of network #2

		// Set interconnections between the two networks.
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[0]);
		DecrementalMCCAlgorithm::interconnect(nodes1[1], nodes2[1]);

		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);

		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[1], edges1[0]);
		SpanningForest::connect(nodes2[1], nodes2[2], edges1[1]);

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1));
		REQUIRE(cluster(nodes1[2]) == nodeSet(nodes1, 2));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1));
		REQUIRE(cluster(nodes2[2]) == nodeSet(nodes2, 2));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));
	}

	{
		NodeVector nodes1(3); // nodes of network #1
		NodeVector nodes2(3); // nodes of network #2
		EdgeVector edges1(3); // edges of network #1
		EdgeVector edges2(3); // edges of network #2

		// Set interconnections between the two networks.
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[0]);
		DecrementalMCCAlgorithm::interconnect(nodes1[1], nodes2[1]);
		DecrementalMCCAlgorithm::interconnect(nodes1[1], nodes2[2]);

		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);

		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[1], edges1[0]);
		SpanningForest::connect(nodes2[1], nodes2[2], edges1[1]);

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1));
		REQUIRE(cluster(nodes1[2]) == nodeSet(nodes1, 2));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1, 2));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));
	}

	{
		NodeVector nodes1(2); // nodes of network #1
		NodeVector nodes2(3); // nodes of network #2
		EdgeVector edges1(3); // edges of network #1
		EdgeVector edges2(3); // edges of network #2

		// Set interconnections between the two networks.
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[0]);
		DecrementalMCCAlgorithm::interconnect(nodes1[1], nodes2[1]);
		DecrementalMCCAlgorithm::interconnect(nodes1[1], nodes2[2]);

		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);

		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[1], edges1[0]);
		SpanningForest::connect(nodes2[1], nodes2[2], edges1[1]);

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1, 2));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));
	}

	{
		NodeVector nodes1(4); // nodes of network #1
		NodeVector nodes2(4); // nodes of network #2
		EdgeVector edges1(6); // edges of network #1
		EdgeVector edges2(6); // edges of network #2

		// Set interconnections between the two networks.
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[0]);
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[1]);
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[2]);
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[3]);

		// Configuration of network #1.
		SpanningForest::connect(nodes1[1], nodes1[2], edges1[0]);
		SpanningForest::connect(nodes1[2], nodes1[3], edges1[1]);
		SpanningForest::connect(nodes1[3], nodes1[1], edges1[2]);

		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[1], edges2[0]);
		SpanningForest::connect(nodes2[1], nodes2[2], edges2[1]);
		SpanningForest::connect(nodes2[2], nodes2[3], edges2[2]);

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0));
		REQUIRE(cluster(nodes1[1]) == nodeSet(nodes1, 1));
		REQUIRE(cluster(nodes1[2]) == nodeSet(nodes1, 2));
		REQUIRE(cluster(nodes1[3]) == nodeSet(nodes1, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));

		// Tracing the distribution of MCC sizes by deleting edges one by one.
		DecrementalMCCAlgorithm::disconnect(edges1[0], reprs1, reprs2);
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));

		DecrementalMCCAlgorithm::disconnect(edges1[1], reprs1, reprs2);
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));

		DecrementalMCCAlgorithm::disconnect(edges1[2], reprs1, reprs2);
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));

		DecrementalMCCAlgorithm::disconnect(edges2[2], reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 0));
		REQUIRE(checkReprSet(nodes2, reprs2, 0));
	}

	{
		NodeVector nodes1(1); // nodes of network #1
		NodeVector nodes2(4); // nodes of network #2
		EdgeVector edges1(0); // edges of network #1
		EdgeVector edges2(6); // edges of network #2

		// Set interconnections between the two networks.
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[0]);
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[1]);
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[2]);
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[3]);

		// Configuration of network #1.

		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[1], edges2[0]);
		SpanningForest::connect(nodes2[1], nodes2[2], edges2[1]);
		SpanningForest::connect(nodes2[2], nodes2[3], edges2[2]);

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));

		// Tracing the distribution of MCC sizes by deleting edges one by one.
		DecrementalMCCAlgorithm::disconnect(edges2[2], reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 0));
		REQUIRE(checkReprSet(nodes2, reprs2, 0));
	}

	{
		NodeVector nodes1(4); // nodes of network #1
		NodeVector nodes2(4); // nodes of network #2
		EdgeVector edges1(6); // edges of network #1
		EdgeVector edges2(6); // edges of network #2

		// Set interconnections between the two networks.
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[0]);
		DecrementalMCCAlgorithm::interconnect(nodes1[0], nodes2[1]);
		DecrementalMCCAlgorithm::interconnect(nodes1[1], nodes2[2]);
		DecrementalMCCAlgorithm::interconnect(nodes1[1], nodes2[3]);

		// Configuration of network #1.
		SpanningForest::connect(nodes1[0], nodes1[1], edges1[0]);
		SpanningForest::connect(nodes1[2], nodes1[3], edges1[1]);

		// Configuration of network #2.
		SpanningForest::connect(nodes2[0], nodes2[1], edges2[0]);
		SpanningForest::connect(nodes2[1], nodes2[2], edges2[1]);
		SpanningForest::connect(nodes2[2], nodes2[3], edges2[2]);

		// Initial distribution of MCC sizes.
		ClusterRepresentativeSet reprs1, reprs2;
		DecrementalMCCAlgorithm::initialize(nodes1, nodes2, reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0, 1));
		REQUIRE(cluster(nodes1[2]) == nodeSet(nodes1, 2));
		REQUIRE(cluster(nodes1[3]) == nodeSet(nodes1, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));

		// Tracing the distribution of MCC sizes by deleting edges one by one.
		DecrementalMCCAlgorithm::disconnect(edges1[0], reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0));
		REQUIRE(cluster(nodes1[1]) == nodeSet(nodes1, 1));
		REQUIRE(cluster(nodes1[2]) == nodeSet(nodes1, 2));
		REQUIRE(cluster(nodes1[3]) == nodeSet(nodes1, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0, 1));
		REQUIRE(cluster(nodes2[2]) == nodeSet(nodes2, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 2));
		REQUIRE(checkReprSet(nodes2, reprs2, 2));

		DecrementalMCCAlgorithm::disconnect(edges2[0], reprs1, reprs2);
		REQUIRE(cluster(nodes1[0]) == nodeSet(nodes1, 0));
		REQUIRE(cluster(nodes1[1]) == nodeSet(nodes1, 1));
		REQUIRE(cluster(nodes1[2]) == nodeSet(nodes1, 2));
		REQUIRE(cluster(nodes1[3]) == nodeSet(nodes1, 3));
		REQUIRE(cluster(nodes2[0]) == nodeSet(nodes2, 0));
		REQUIRE(cluster(nodes2[1]) == nodeSet(nodes2, 1));
		REQUIRE(cluster(nodes2[2]) == nodeSet(nodes2, 2, 3));
		REQUIRE(checkReprSet(nodes1, reprs1, 1));
		REQUIRE(checkReprSet(nodes2, reprs2, 1));

		DecrementalMCCAlgorithm::disconnect(edges2[2], reprs1, reprs2);
		REQUIRE(ensureTotallyDisconnected(nodes1, nodes2));
		REQUIRE(checkReprSet(nodes1, reprs1, 0));
		REQUIRE(checkReprSet(nodes2, reprs2, 0));
	}
}