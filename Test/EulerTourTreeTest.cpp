#define CATCH_CONFIG_MAIN
#include "catch/catch.hpp"
#include <memory>
#include "EulerTourTree.hpp"

using namespace std;
using namespace Snu::Cnrc::EulerTourTree;

struct MyNode;
struct MyEdge;

struct MyNode : public EulerTourTreeNodeMixin<MyNode, MyEdge> {
	
	//for convenience of test
	bool operator==(const MyNode& other) const {
		return this == &other;
	}
	
	//for convenience of test
	bool operator!=(const MyNode& other) const {
		return this != &other;
	}
	
	//for convenience of test.
	bool operator<(const MyNode& other) const {
		return this < &other;
	}
	
	//for convenience of test.
	unsigned int name;
};

struct MyEdge : public EulerTourTreeEdgeMixin<MyNode, MyEdge> {
	
	//for convenience of test.
	unsigned int node1;
	unsigned int node2;
};

using ETT = EulerTourTreeAlgorithm<MyNode, MyEdge>;

TEST_CASE("TestSingleNodeTree") {
	MyNode node;
	node.name = 1234;
	REQUIRE(ETT::hasPath(node, node));
	REQUIRE(ETT::isClusterRep(node));
	REQUIRE(ETT::clusterSize(node) == 1);
	REQUIRE(ETT::findClusterRep(node) == node);
	
	{
		ETT::NodeContainerView seq = ETT::nodeContainerView(node);
		ETT::NodeContainerView::iterator itr = seq.begin();
		REQUIRE(*itr == node);
		REQUIRE(itr->name == 1234);
		REQUIRE(ETT::isClusterRep(*itr));
		itr->name = 4321;
		REQUIRE(*itr++ == node);
		REQUIRE(itr == seq.end());
	}

	{
		ETT::ConstNodeContainerView seq = ETT::nodeContainerView(static_cast<const MyNode&>(node));
		ETT::ConstNodeContainerView::const_iterator itr = seq.begin();
		REQUIRE(*itr == node);
		REQUIRE(itr->name == 4321);
		REQUIRE(ETT::isClusterRep(*itr));
		REQUIRE(*itr++ == node);
		REQUIRE(itr == seq.end());
	}
}

set<const MyNode*> toNodeSet(const MyNode& n) {
	set<const MyNode*> s;
	for(const MyNode& m : ETT::nodeContainerView(n)) {
		REQUIRE(s.insert(&m).second);
	}
	return s;
}

set<const MyEdge*> toEdgeSet(const MyNode& n) {
	set<const MyEdge*> s;
	for(const MyEdge& e : ETT::edgeContainerView(n)) {
		REQUIRE(s.insert(&e).second);
	}
	return s;
}

MyNode* makeNodes(unsigned int n) {
	MyNode* const nodes = new MyNode[n];
	for(unsigned int i = 0;  i < n; ++i) {
		nodes[i].name = i;
	}
	return nodes;
}

MyEdge* makeEdges(unsigned int n) {
	MyEdge* const edges = new MyEdge[n];
	for(unsigned int i = 0;  i < n; ++i) {
		edges[i].node1 = i;
		edges[i].node2 = i + 1;
	}
	return edges;
}

void connectRange(
	MyNode* const nodes,
	MyEdge* const edges,
	const unsigned int from, 
	const unsigned to
) {
	for(unsigned int i = from;  i < to - 1; ++i) {
		ETT::connect(nodes[i], nodes[i + 1], edges[i]);
	}
}

void assertRangeConnected(
	const MyNode* const nodes,
	const MyEdge* const edges, 
	const unsigned int from, 
	const unsigned to
) {
	const MyNode& r = *find_if(nodes + from, nodes + to, [](const MyNode& n) {
		return ETT::isClusterRep(n);
	});
	set<const MyNode*> ns;
	set<const MyEdge*> es;
	for(unsigned int i = from; i < to; ++i) {
		ns.insert(nodes + i);
	}
	for(unsigned int i = from; i < to - 1; ++i) {
		es.insert(edges + i);
	}
	for(unsigned int i = from; i < to; ++i) {
		if(nodes[i] != r) {
			REQUIRE(! ETT::isClusterRep(nodes[i]));
		}
		REQUIRE(ETT::findClusterRep(nodes[i]) == r);
		REQUIRE(ETT::clusterSize(nodes[i]) == to - from);
		for(unsigned int j = from; j < to; ++j) {
			REQUIRE(ETT::hasPath(nodes[i], nodes[j]));
		}
		REQUIRE(toNodeSet(nodes[i]) == ns);
		REQUIRE(toEdgeSet(nodes[i]) == es);
	}
}

void assertRangeNotConnected(
	const MyNode* const nodes, 
	const unsigned int from1, 
	const unsigned int to1,
	const unsigned int from2, 
	const unsigned int to2
) {
	for(unsigned int i = from1; i < to1; ++i) {
		for(unsigned int j = from2; j < to2; ++j) {
			REQUIRE(! ETT::hasPath(nodes[i], nodes[j]));
		}
	}
}

TEST_CASE("TestTwoNodesTree") {
	MyNode nodes[2];
	MyEdge edge;
	assertRangeNotConnected(nodes, 0, 1, 1, 2);
	ETT::connect(nodes[0], nodes[1], edge);
	assertRangeConnected(nodes, &edge, 0, 2);
	ETT::disconnect(edge);
	assertRangeNotConnected(nodes, 0, 1, 1, 2);
	assertRangeConnected(nodes, &edge, 0, 1);
	assertRangeConnected(nodes, &edge, 1, 2);
}

TEST_CASE("TestNameConfliction") {

	struct MyNode;
	struct MyEdge;
	
	struct MyNode : public EulerTourTreeNodeMixin<MyNode, MyEdge> {
		
		bool operator==(MyNode& other) {
			return this == &other;
		}
		
		int activeOccurrence;
	};
	
	struct MyEdge : public EulerTourTreeEdgeMixin<MyNode, MyEdge> {
		
		bool operator==(MyEdge& other) {
			return this == &other;
		}
		
		int occurrence1;
		int occurrence2;
		int occurrence3;
		int occurrence4;
		
	};
	
	using ETT = EulerTourTreeAlgorithm<MyNode, MyEdge>;

	MyNode nodes[2];
	REQUIRE(*ETT::nodeContainerView(nodes[0]).begin() == nodes[0]);
	REQUIRE(*ETT::nodeContainerView(const_cast<MyNode&>(nodes[0])).begin() == nodes[0]);
	MyEdge edge;
	ETT::connect(nodes[0], nodes[1], edge);
	REQUIRE(ETT::hasPath(nodes[0], nodes[1]));
	REQUIRE(ETT::clusterSize(nodes[0]) == 2);
	REQUIRE(ETT::clusterSize(nodes[1]) == 2);
	REQUIRE(*ETT::edgeContainerView(nodes[0]).begin() == edge);
	REQUIRE(*ETT::edgeContainerView(const_cast<MyNode&>(nodes[0])).begin() == edge);
	ETT::disconnect(edge);
	REQUIRE(! ETT::hasPath(nodes[0], nodes[1]));
	REQUIRE(ETT::clusterSize(nodes[0]) == 1);
	REQUIRE(ETT::isClusterRep(nodes[0]));
	REQUIRE(&ETT::findClusterRep(nodes[0]) == nodes + 0);
	REQUIRE(ETT::clusterSize(nodes[1]) == 1);
	REQUIRE(ETT::isClusterRep(static_cast<const MyNode&>(nodes[1])));
	REQUIRE(&ETT::findClusterRep(static_cast<const MyNode&>(nodes[1])) == nodes + 1);
}

TEST_CASE("TestManyNodesTreeA") {
	for(unsigned int n = 3; n < 16; ++n) {
		unique_ptr<MyNode[]> nodes(makeNodes(n));
		unique_ptr<MyEdge[]> edges(makeEdges(n - 1));
		for(unsigned int i = 0; i < n - 1; ++i) {
			ETT::connect(nodes[i], nodes[i + 1], edges[i]);
			assertRangeConnected(nodes.get(), edges.get(), 0, i + 2);
			assertRangeNotConnected(nodes.get(), 0, i + 2, i + 2, n);
		}
		for(unsigned int i = n - 1; i > 0; --i) {
			ETT::disconnect(edges[i - 1]); //deleting the edge between i - 1 and i.
			assertRangeConnected(nodes.get(), edges.get(), 0, i);
			REQUIRE(toNodeSet(nodes[i]) == set<const MyNode*>{nodes.get() + i});
			assertRangeNotConnected(nodes.get(), 0, i, i, n);
			REQUIRE(ETT::isClusterRep(nodes[i]));
			REQUIRE(ETT::clusterSize(nodes[i]) == 1);
		}
	}
}

TEST_CASE("TestManyNodesTreeB") {
	const unsigned int n = 32;
	for(unsigned int i = 1; i < n; ++i) {
		unique_ptr<MyNode[]> nodes(makeNodes(n));
		unique_ptr<MyEdge[]> edges(makeEdges(n - 1));
		connectRange(nodes.get(), edges.get(), 0, i);
		connectRange(nodes.get(), edges.get(), i, n);
		assertRangeConnected(nodes.get(), edges.get(), 0, i);
		assertRangeConnected(nodes.get(), edges.get(), i, n);
		assertRangeNotConnected(nodes.get(), 0, i, i, n);
		ETT::connect(nodes[i - 1], nodes[i], edges[i - 1]);
		assertRangeConnected(nodes.get(), edges.get(), 0, n);
		for(unsigned int j = 0; j < n - 1; ++j) {
			ETT::disconnect(edges[j]); //deleting the edge between j and j + 1.
			assertRangeConnected(nodes.get(), edges.get(), 0, j + 1);
			assertRangeConnected(nodes.get(), edges.get(), j + 1, n);
			assertRangeNotConnected(nodes.get(), 0, j + 1, j + 1, n);
			ETT::connect(nodes[j], nodes[j + 1], edges[j]);
			assertRangeConnected(nodes.get(), edges.get(), 0, n);
		}
	}
}
