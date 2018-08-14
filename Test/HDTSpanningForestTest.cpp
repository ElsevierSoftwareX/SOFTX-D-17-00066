#define CATCH_CONFIG_MAIN
#include "catch/catch.hpp"
#include <memory>
#include <set>
#include <map>
#include <utility>
#include "HDTSpanningForest.hpp"

using namespace std;
using namespace Snu::Cnrc::HDTSpanningForest;

struct MyNode;
struct MyEdge;

struct MyNode : public HDTSpanningForestNodeMixin<MyNode, MyEdge> {
	
	MyNode()
	: HDTSpanningForestNodeMixin() {
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

struct MyEdge : public HDTSpanningForestEdgeMixin<MyNode, MyEdge> {
};

using HDT = HDTSpanningForestAlgorithm<MyNode, MyEdge>;

TEST_CASE("TestSingleNodeGraph") {
	MyNode n;
	n.name = 1234;
	REQUIRE(HDT::hasPath(n, n));
	REQUIRE(HDT::isClusterRep(n));
	REQUIRE(HDT::findClusterRep(n) == n);
	
	{
		HDT::Cluster c = HDT::cluster(n);
		REQUIRE(c.size() == 1);
		HDT::Cluster::iterator itr = c.begin();
		REQUIRE(*itr == n);
		REQUIRE(itr->name == 1234);
		itr->name = 4321;
		REQUIRE(*itr++ == n);
		REQUIRE(itr == c.end());
	}
	
	{
		HDT::ConstCluster c = HDT::cluster(static_cast<const MyNode&>(n));
		REQUIRE(c.size() == 1);
		HDT::Cluster::const_iterator itr = c.begin();
		REQUIRE(*itr == n);
		REQUIRE(itr->name == 4321);
		REQUIRE(*itr++ == n);
		REQUIRE(itr == c.end());
	}

	{
		REQUIRE(HDT::edges(n).begin() == HDT::edges(n).end());
		REQUIRE(HDT::edges(n).cbegin() == HDT::edges(n).cend());
		REQUIRE(HDT::edges(static_cast<const MyNode&>(n)).begin() == 
					  HDT::edges(static_cast<const MyNode&>(n)).end()
		);
	}
}

set<MyNode*> cluster(MyNode& n) {
	set<MyNode*> s;
	for(MyNode& m : HDT::cluster(n)) {
		REQUIRE(s.insert(&m).second);
	}
	return s;
}

set<const MyNode*> cluster(const MyNode& n) {
	set<const MyNode*> s;
	for(const MyNode& m : HDT::cluster(n)) {
		REQUIRE(s.insert(&m).second);
	}
	return s;
}

set<const MyNode*> constCluster(MyNode& n) {
	set<const MyNode*> s;
	for(const MyNode& m : HDT::constCluster(n)) {
		REQUIRE(s.insert(&m).second);
	}
	return s;
}

set<const MyEdge*> edges(const MyNode& n) {
	set<const MyEdge*> s;
	for(const MyEdge& m : HDT::edges(n)) {
		REQUIRE(s.insert(&m).second);
	}
	return s;
}

set<MyEdge*> edges(MyNode& n) {
	set<MyEdge*> s;
	for(MyEdge& m : HDT::edges(n)) {
		REQUIRE(s.insert(&m).second);
	}
	return s;
}

set<const MyEdge*> constEdges(MyNode& n) {
	set<const MyEdge*> s;
	for(const MyEdge& m : HDT::constEdges(n)) {
		REQUIRE(s.insert(&m).second);
	}
	return s;
}

TEST_CASE("TestTwoNodesGraph") {
	MyNode nodes[2];
	MyEdge edge;
	REQUIRE(! HDT::isValid(edge));
	REQUIRE(! HDT::hasPath(nodes[0], nodes[1]));
	REQUIRE(HDT::connect(nodes[0], nodes[1], edge));
	REQUIRE(HDT::isValid(edge));
	REQUIRE(HDT::hasPath(nodes[0], nodes[1]));
	REQUIRE(HDT::clusterSize(nodes[0]) == 2);
	REQUIRE(HDT::clusterSize(nodes[1]) == 2);
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
	REQUIRE(HDT::disconnect(edge));
	REQUIRE(! HDT::isValid(edge));
	REQUIRE(! HDT::hasPath(nodes[0], nodes[1]));
	REQUIRE(cluster(nodes[0]) == set<MyNode*>({nodes}));
	REQUIRE(cluster(nodes[1]) == set<MyNode*>({nodes + 1}));
	REQUIRE(cluster(static_cast<const MyNode&>(nodes[0])) == set<const MyNode*>({nodes}));
	REQUIRE(cluster(static_cast<const MyNode&>(nodes[1])) == set<const MyNode*>({nodes + 1}));
	REQUIRE(edges(nodes[0]) == set<MyEdge*>());
	REQUIRE(edges(nodes[1]) == set<MyEdge*>());
	REQUIRE(edges(static_cast<const MyNode&>(nodes[0])) == set<const MyEdge*>());
	REQUIRE(edges(static_cast<const MyNode&>(nodes[1])) == set<const MyEdge*>());
	REQUIRE(HDT::clusterSize(nodes[0]) == 1);
	REQUIRE(HDT::clusterSize(nodes[1]) == 1);
}

TEST_CASE("TestNameConfliction") {
	
	struct MyNode;
	struct MyEdge;

	struct MyNode : public HDTSpanningForestNodeMixin<MyNode, MyEdge> {
		int levelNodes;
	};
	
	struct MyEdge : public HDTSpanningForestEdgeMixin<MyNode, MyEdge> {
		double node1;
		double node2;
		double level;
		double isTreeEdge;
		double levelEdges;
	};
	
	using HDT = HDTSpanningForestAlgorithm<MyNode, MyEdge>;

	MyNode nodes[2];
	MyEdge edge;
	
	REQUIRE(! HDT::hasPath(nodes[0], nodes[1]));
	REQUIRE(HDT::connect(nodes[0], nodes[1], edge));
	REQUIRE(HDT::hasPath(nodes[0], nodes[1]));
	REQUIRE(HDT::cluster(nodes[0]).size() == 2);
	REQUIRE(HDT::clusterSize(nodes[1]) == 2);
	REQUIRE(HDT::isClusterRep(HDT::findClusterRep(nodes[0])));
	REQUIRE(HDT::disconnect(edge));
	REQUIRE(! HDT::hasPath(nodes[0], nodes[1]));
	REQUIRE(HDT::cluster(static_cast<const MyNode&>(nodes[0])).size() == 1);
	REQUIRE(HDT::cluster(static_cast<const MyNode&>(nodes[1])).size() == 1);
	REQUIRE(HDT::isClusterRep(HDT::findClusterRep(static_cast<const MyNode&>(nodes[0]))));
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
		return HDT::isClusterRep(n);
	});
	set<const MyNode*> ss;
	for(unsigned int i = from; i < to; ++i) {
		ss.insert(nodes + i);
	}
	for(unsigned int i = from; i < to; ++i) {
		if(nodes[i] != r) {
			REQUIRE(! HDT::isClusterRep(nodes[i]));
		
		REQUIRE(HDT::findClusterRep(nodes[i]) == r);}
		for(unsigned int j = from; j < to; ++j) {
			REQUIRE(HDT::hasPath(nodes[i], nodes[j]));
		}
		REQUIRE(HDT::cluster(nodes[i]).size() == to - from);
		REQUIRE(HDT::clusterSize(nodes[i]) == to - from);
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
			REQUIRE(! HDT::hasPath(nodes[i], nodes[j]));
		}
	}
}

TEST_CASE("TestManyNodesGraph1") {
	
	using NodeNamePair = pair<unsigned int, unsigned int>;
	
	for(unsigned int n = 3; n <= 32; ++n) {
		unique_ptr<MyNode[]> nodes(makeNodes(n));
		map<NodeNamePair, unique_ptr<MyEdge>> edges;
		for(unsigned int i = 0; i < n; ++i) {
			for(unsigned int j = i / 2; j < i; ++j) {
				edges[NodeNamePair(i, j)] = unique_ptr<MyEdge>(new MyEdge());
				bool merged = HDT::connect(nodes[i], nodes[j], *edges[NodeNamePair(i, j)]);
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
					REQUIRE(! HDT::disconnect(*e));
					REQUIRE(! HDT::isValid(*e));
					edges.erase(NodeNamePair(i, j));
					assertRangeConnected(nodes.get(), 0, n);
				}
			}
		}
		for(unsigned int i = 0; i < n - 1; ++i) {
			REQUIRE(HDT::disconnect(*edges[NodeNamePair(i + 1, i)]));
			REQUIRE(HDT::cluster(nodes[i]).size() == 1);
			REQUIRE(HDT::clusterSize(nodes[i]) == 1);
			REQUIRE(cluster(nodes[i]) == set<MyNode*>{nodes.get() + i});
			assertRangeConnected(nodes.get(), i + 1, n);
			REQUIRE(HDT::cluster(nodes[i + 1]).size() == n - i - 1);
			REQUIRE(HDT::clusterSize(nodes[i + 1]) == n - i - 1);
		}
	}
}

TEST_CASE("TestManyNodesGraph2") {
	
	using NodeNamePair = pair<unsigned int, unsigned int>;
	
	for(unsigned int n = 3; n <= 32; ++n) {
		unique_ptr<MyNode[]> nodes(makeNodes(n));
		map<NodeNamePair, unique_ptr<MyEdge>> edges;
		for(unsigned int i = 0; i < n; ++i) {
			for(unsigned int j = i / 2; j < i; ++j) {
				edges[NodeNamePair(i, j)] = unique_ptr<MyEdge>(new MyEdge());
				bool merged = HDT::connect(nodes[i], nodes[j], *edges[NodeNamePair(i, j)]);
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
					REQUIRE(! HDT::disconnect(*e));
					REQUIRE(! HDT::isValid(*e));
					edges.erase(NodeNamePair(i, j));
					assertRangeConnected(nodes.get(), 0, n);
				}
			}
		}
		for(unsigned int i = 0; i < n - 1; ++i) {
			REQUIRE(HDT::disconnect(*edges[NodeNamePair(i + 1, i)]));
			REQUIRE(HDT::cluster(nodes[i]).size() == 1);
			REQUIRE(HDT::clusterSize(nodes[i]) == 1);
			REQUIRE(cluster(nodes[i]) == set<MyNode*>{nodes.get() + i});
			assertRangeConnected(nodes.get(), i + 1, n);
			REQUIRE(HDT::cluster(nodes[i + 1]).size() == n - i - 1);
			REQUIRE(HDT::clusterSize(nodes[i + 1]) == n - i - 1);
		}
	}
}
