#define CATCH_CONFIG_MAIN
#include "catch/catch.hpp"
#include <memory>
#include <set>
#include <map>
#include <utility>
#include "EulerTourTreeSpanningForest.hpp"

using namespace std;
using namespace Snu::Cnrc::EulerTourTreeSpanningForest;

struct MyNode;
struct MyEdge;

struct MyNode : public EulerTourTreeSpanningForestNodeMixin<MyNode, MyEdge> {
	
	MyNode()
	: EulerTourTreeSpanningForestNodeMixin() {
	}
	
	//to test with a minimal restricted class
	MyNode(const MyNode&) = delete;
	
	//to test with a minimal restricted class
	MyNode(MyNode&&) = delete;
	
	//to test with a minimal restricted class
	MyNode& operator=(const MyNode&) = delete;
	
	//to test with a minimal restricted class
	MyNode& operator=(MyNode&&) = delete;
	
	//for convenience of test
	bool operator==(const MyNode& other) const {
		return this == &other;
	}
	
	//for convenience of test
	bool operator!=(const MyNode& other) const {
		return this != &other;
	}
	
	unsigned int name;
};

struct MyEdge : public EulerTourTreeSpanningForestEdgeMixin<MyNode, MyEdge> {
};

using ETTSpanningForest = EulerTourTreeSpanningForestAlgorithm<MyNode, MyEdge>;

TEST_CASE("TestSingleNodeGraph") {
	MyNode n;
	n.name = 1234;
	REQUIRE(ETTSpanningForest::hasPath(n, n));
	REQUIRE(ETTSpanningForest::isClusterRep(n));
	REQUIRE(ETTSpanningForest::findClusterRep(n) == n);
	
	{
		ETTSpanningForest::Cluster c = ETTSpanningForest::cluster(n);
		REQUIRE(c.size() == 1);
		ETTSpanningForest::Cluster::iterator itr = c.begin();
		REQUIRE(*itr == n);
		REQUIRE(itr->name == 1234);
		itr->name = 4321;
		REQUIRE(*itr++ == n);
		REQUIRE(itr == c.end());
	}
	
	{
		ETTSpanningForest::ConstCluster c = ETTSpanningForest::cluster(static_cast<const MyNode&>(n));
		REQUIRE(c.size() == 1);
		ETTSpanningForest::Cluster::const_iterator itr = c.begin();
		REQUIRE(*itr == n);
		REQUIRE(itr->name == 4321);
		REQUIRE(*itr++ == n);
		REQUIRE(itr == c.end());
	}

	{
		REQUIRE(ETTSpanningForest::edges(n).begin() == ETTSpanningForest::edges(n).end());
		REQUIRE(ETTSpanningForest::edges(n).cbegin() == ETTSpanningForest::edges(n).cend());
		REQUIRE(ETTSpanningForest::edges(static_cast<const MyNode&>(n)).begin() == 
					  ETTSpanningForest::edges(static_cast<const MyNode&>(n)).end()
		);
	}
}

set<const MyNode*> cluster(const MyNode& n) {
	set<const MyNode*> s;
	for(const MyNode& m : ETTSpanningForest::cluster(n)) {
		REQUIRE(s.insert(&m).second);
	}
	return s;
}

set<MyNode*> cluster(MyNode& n) {
	set<MyNode*> s;
	for(MyNode& m : ETTSpanningForest::cluster(n)) {
		REQUIRE(s.insert(&m).second);
	}
	return s;
}

set<const MyNode*> constCluster(MyNode& n) {
	set<const MyNode*> s;
	for(const MyNode& m : ETTSpanningForest::constCluster(n)) {
		REQUIRE(s.insert(&m).second);
	}
	return s;
}

set<const MyEdge*> edges(const MyNode& n) {
	set<const MyEdge*> s;
	for(const MyEdge& m : ETTSpanningForest::edges(n)) {
		REQUIRE(s.insert(&m).second);
	}
	return s;
}

set<MyEdge*> edges(MyNode& n) {
	set<MyEdge*> s;
	for(MyEdge& m : ETTSpanningForest::edges(n)) {
		REQUIRE(s.insert(&m).second);
	}
	return s;
}

set<const MyEdge*> constEdges(MyNode& n) {
	set<const MyEdge*> s;
	for(const MyEdge& m : ETTSpanningForest::constEdges(n)) {
		REQUIRE(s.insert(&m).second);
	}
	return s;
}

TEST_CASE("TestTwoNodesGraph") {
	MyNode nodes[2];
	MyEdge edge;
	REQUIRE(! ETTSpanningForest::isValid(edge));
	REQUIRE(! ETTSpanningForest::hasPath(nodes[0], nodes[1]));
	REQUIRE(ETTSpanningForest::connect(nodes[0], nodes[1], edge));
	REQUIRE(ETTSpanningForest::isValid(edge));
	REQUIRE(ETTSpanningForest::hasPath(nodes[0], nodes[1]));
	REQUIRE(ETTSpanningForest::clusterSize(nodes[0]) == 2);
	REQUIRE(ETTSpanningForest::clusterSize(nodes[1]) == 2);
	REQUIRE(cluster(nodes[0]) == set<MyNode*>({nodes, nodes + 1}));
	REQUIRE(cluster(nodes[1]) == set<MyNode*>({nodes, nodes + 1}));
	REQUIRE(cluster(static_cast<const MyNode&>(nodes[0])) == set<const MyNode*>({nodes, nodes + 1}));
	REQUIRE(cluster(static_cast<const MyNode&>(nodes[1])) == set<const MyNode*>({nodes, nodes + 1}));
	REQUIRE(constCluster(nodes[0]) == set<const MyNode*>({nodes, nodes + 1}));
	REQUIRE(constCluster(nodes[1]) == set<const MyNode*>({nodes, nodes + 1}));
	REQUIRE(edges(nodes[0]) == set<MyEdge*>({&edge}));
	REQUIRE(edges(nodes[1]) == set<MyEdge*>({&edge}));
	REQUIRE(edges(static_cast<const MyNode&>(nodes[0])) == set<const MyEdge*>({&edge}));
	REQUIRE(edges(static_cast<const MyNode&>(nodes[1])) == set<const MyEdge*>({&edge}));
	REQUIRE(constEdges(nodes[0]) == set<const MyEdge*>({&edge}));
	REQUIRE(constEdges(nodes[1]) == set<const MyEdge*>({&edge}));
	REQUIRE(ETTSpanningForest::disconnect(edge));
	REQUIRE(! ETTSpanningForest::isValid(edge));
	REQUIRE(! ETTSpanningForest::hasPath(nodes[0], nodes[1]));
	REQUIRE(cluster(nodes[0]) == set<MyNode*>({nodes}));
	REQUIRE(cluster(nodes[1]) == set<MyNode*>({nodes + 1}));
	REQUIRE(cluster(static_cast<const MyNode&>(nodes[0])) == set<const MyNode*>({nodes}));
	REQUIRE(cluster(static_cast<const MyNode&>(nodes[1])) == set<const MyNode*>({nodes + 1}));
	REQUIRE(edges(nodes[0]) == set<MyEdge*>());
	REQUIRE(edges(nodes[1]) == set<MyEdge*>());
	REQUIRE(edges(static_cast<const MyNode&>(nodes[0])) == set<const MyEdge*>());
	REQUIRE(edges(static_cast<const MyNode&>(nodes[1])) == set<const MyEdge*>());
	REQUIRE(ETTSpanningForest::clusterSize(nodes[0]) == 1);
	REQUIRE(ETTSpanningForest::clusterSize(nodes[1]) == 1);
}

TEST_CASE("Test Name Confliction") {
	
	struct MyNode;
	struct MyEdge;

	struct MyNode : public EulerTourTreeSpanningForestNodeMixin<MyNode, MyEdge> {
		int edges;
	};
	
	struct MyEdge : public EulerTourTreeSpanningForestEdgeMixin<MyNode, MyEdge> {
		double node1;
		double node2;
		double isTreeEdge;
	};
	
	using ETTSpanningForest = EulerTourTreeSpanningForestAlgorithm<MyNode, MyEdge>;

	MyNode nodes[2];
	MyEdge edge;
	
	REQUIRE(! ETTSpanningForest::hasPath(nodes[0], nodes[1]));
	REQUIRE(ETTSpanningForest::connect(nodes[0], nodes[1], edge));
	REQUIRE(ETTSpanningForest::hasPath(nodes[0], nodes[1]));
	REQUIRE(ETTSpanningForest::cluster(nodes[0]).size() == 2);
	REQUIRE(ETTSpanningForest::clusterSize(nodes[1]) == 2);
	REQUIRE(ETTSpanningForest::isClusterRep(ETTSpanningForest::findClusterRep(nodes[0])));
	REQUIRE(ETTSpanningForest::disconnect(edge));
	REQUIRE(! ETTSpanningForest::hasPath(nodes[0], nodes[1]));
	REQUIRE(ETTSpanningForest::cluster(nodes[1]).size() == 1);
	REQUIRE(ETTSpanningForest::isClusterRep(ETTSpanningForest::findClusterRep(nodes[0])));
}

MyNode* makeNodes(unsigned int n) {
	MyNode* const nodes = new MyNode[n];
	for(unsigned int i = 0;  i < n; ++i) {
		nodes[i].name = i;
	}
	return nodes;
}

void assertRangeConnected(
	const MyNode* const nodes, 
	const unsigned int from, 
	const unsigned to
) {
	const MyNode& r = *find_if(nodes + from, nodes + to, [](const MyNode& n) {
		return ETTSpanningForest::isClusterRep(n);
	});
	set<const MyNode*> ss;
	for(unsigned int i = from; i < to; ++i) {
		ss.insert(nodes + i);
	}
	for(unsigned int i = from; i < to; ++i) {
		if(nodes[i] != r) {
			REQUIRE(! ETTSpanningForest::isClusterRep(nodes[i]));
		
		REQUIRE(ETTSpanningForest::findClusterRep(nodes[i]) == r);}
		for(unsigned int j = from; j < to; ++j) {
			REQUIRE(ETTSpanningForest::hasPath(nodes[i], nodes[j]));
		}
		REQUIRE(ETTSpanningForest::cluster(nodes[i]).size() == to - from);
		REQUIRE(ETTSpanningForest::clusterSize(nodes[i]) == to - from);
		REQUIRE(cluster(nodes[i]) == ss);
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
			REQUIRE(! ETTSpanningForest::hasPath(nodes[i], nodes[j]));
		}
	}
}

TEST_CASE("TestManyNodesGraph1") {
	
	using NodeNamePair = pair<unsigned int, unsigned int>;
	
	for(unsigned int n = 3; n < 32; ++n) {
		unique_ptr<MyNode[]> nodes(makeNodes(n));
		map<NodeNamePair, unique_ptr<MyEdge>> edges;
		for(unsigned int i = 0; i < n; ++i) {
			for(unsigned int j = i / 2; j < i; ++j) {
				edges[NodeNamePair(i, j)] = unique_ptr<MyEdge>(new MyEdge());
				bool merged = ETTSpanningForest::connect(nodes[i], nodes[j], *edges[NodeNamePair(i, j)]);
				if(j == i / 2) {
					REQUIRE(merged);
					assertRangeConnected(nodes.get(), 0, i + 1);
				} else {
					REQUIRE(! merged);
				}
			}
		}
		assertRangeConnected(nodes.get(), 0, n);
		for(unsigned int i = 0; i < n; ++i) {
			for(unsigned int j = i / 2; j < i; ++j) {
				if(i != j + 1) {
					MyEdge* e = edges[NodeNamePair(i, j)].get();
					REQUIRE(! ETTSpanningForest::disconnect(*e));
					REQUIRE(! ETTSpanningForest::isValid(*e));
					edges.erase(NodeNamePair(i, j));
					assertRangeConnected(nodes.get(), 0, n);
				}
			}
		}
		for(unsigned int i = 0; i < n - 1; ++i) {
			REQUIRE(ETTSpanningForest::disconnect(*edges[NodeNamePair(i + 1, i)]));
			REQUIRE(ETTSpanningForest::cluster(nodes[i]).size() == 1);
			REQUIRE(ETTSpanningForest::clusterSize(nodes[i]) == 1);
			REQUIRE(cluster(nodes[i]) == set<MyNode*>{nodes.get() + i});
			assertRangeConnected(nodes.get(), i + 1, n);
			REQUIRE(ETTSpanningForest::cluster(nodes[i + 1]).size() == n - i - 1);
			REQUIRE(ETTSpanningForest::clusterSize(nodes[i + 1]) == n - i - 1);
		}
	}
}

TEST_CASE("TestManyNodesGraph2") {
	
	using NodeNamePair = pair<unsigned int, unsigned int>;
	
	for(unsigned int n = 3; n < 32; ++n) {
		unique_ptr<MyNode[]> nodes(makeNodes(n));
		map<NodeNamePair, unique_ptr<MyEdge>> edges;
		for(unsigned int i = 0; i < n; ++i) {
			for(unsigned int j = i / 2; j < i; ++j) {
				edges[NodeNamePair(i, j)] = unique_ptr<MyEdge>(new MyEdge());
				bool merged = ETTSpanningForest::connect(nodes[i], nodes[j], *edges[NodeNamePair(i, j)]);
				if(j == i / 2) {
					REQUIRE(merged);
					assertRangeConnected(nodes.get(), 0, i + 1);
				} else {
					REQUIRE(! merged);
				}
			}
		}
		assertRangeConnected(nodes.get(), 0, n);
		for(unsigned int i = 0; i < n; ++i) {
			if(i > 2) {
				for(unsigned int j = i - 2; j >= i / 2; --j) {
					MyEdge* e = edges[NodeNamePair(i, j)].get();
					REQUIRE(! ETTSpanningForest::disconnect(*e));
					REQUIRE(! ETTSpanningForest::isValid(*e));
					edges.erase(NodeNamePair(i, j));
					assertRangeConnected(nodes.get(), 0, n);
				}
			}
		}
		for(unsigned int i = 0; i < n - 1; ++i) {
			REQUIRE(ETTSpanningForest::disconnect(*edges[NodeNamePair(i + 1, i)]));
			REQUIRE(ETTSpanningForest::cluster(nodes[i]).size() == 1);
			REQUIRE(ETTSpanningForest::clusterSize(nodes[i]) == 1);
			REQUIRE(cluster(nodes[i]) == set<MyNode*>{nodes.get() + i});
			assertRangeConnected(nodes.get(), i + 1, n);
			REQUIRE(ETTSpanningForest::cluster(nodes[i + 1]).size() == n - i - 1);
			REQUIRE(ETTSpanningForest::clusterSize(nodes[i + 1]) == n - i - 1);
		}
	}
}
