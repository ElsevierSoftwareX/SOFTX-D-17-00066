/**
 * This program takes three arguments [number of nodes], [initial number of edges] and [final number of edges].
 * It creates a Erdos-Renyi random network with [number of nodes] and [initial number of edges]
 * and then prints the mean degree and the size of the largest cluster while 
 * removing the edges one by one in an randomized order until [final number of edges] remains.
 */

#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <memory>
#include <random>
#include "HDTSpanningForest.hpp"
#include "EulerTourTreeSpanningForest.hpp"

using namespace std;
using namespace boost;

#ifdef CNRC_ENABLE_HDT
	using namespace Snu::Cnrc::HDTSpanningForest;
	using Node = HDTSpanningForestNode;
	using Edge = HDTSpanningForestEdge;
	using SpanningForest = HDTSpanningForestAlgorithm<Node, Edge>;
#else
	using namespace Snu::Cnrc::EulerTourTreeSpanningForest;
	using Node = EulerTourTreeSpanningForestNode;
	using Edge = EulerTourTreeSpanningForestEdge;
	using SpanningForest = EulerTourTreeSpanningForestAlgorithm<Node, Edge>;
#endif

using Size = unsigned int;
using NodeName = Size;
using ClusterSizeDist = map<Size, Size>;
using Graph = map<NodeName, set<NodeName>>;

unique_ptr<Graph> makeERGraph(const Size numberOfNodes, const Size numberOfEdges) {
	Graph* g = new Graph();
	auto rnd = bind(
        uniform_int_distribution<unsigned int>(0, numberOfNodes - 1), 
        default_random_engine(random_device()())
    );
	for(unsigned int i = 0; i < numberOfEdges; ++i) {
		NodeName n1 = rnd();
		NodeName n2 = rnd();
		while(n1 == n2 || (*g)[n1].count(n2)) {
			n1 = rnd();
			n2 = rnd();
		}
		(*g)[n1].insert(n2);
		(*g)[n2].insert(n1);
	}
	return unique_ptr<Graph>(g);
}

void initialize(
	const Graph& g,
	Node* const nodes,
	Edge* const edges
) {
	Size i = 0;
	for(auto p : g) {
		NodeName n = p.first;
		for(auto m : p.second) {
			if(n < m) {
				SpanningForest::connect(nodes[n], nodes[m], edges[i]);
				++i;
			}
		}
	}
}

void putClusterSizes(ClusterSizeDist& cdist, Size sz1, Size sz2) {
	cdist[sz1]++;
	cdist[sz2]++;
	ClusterSizeDist::iterator i = cdist.find(sz1 + sz2);
	if(i->second == 1) {
		cdist.erase(i);
	} else {
		i->second--;
	}
}

void disconnect(
	ClusterSizeDist& cdist,
	Edge& e
) {
	Node& n1 = SpanningForest::node1(e);
	Node& n2 = SpanningForest::node2(e);
	if(! SpanningForest::disconnect(e)) {
		return;
	}
	const Size sz1 = SpanningForest::clusterSize(n1);
	const Size sz2 = SpanningForest::clusterSize(n2);
	putClusterSizes(cdist, sz1, sz2);
}

Size largestClusterSize(const Node* const nodes, Size numberOfNodes) {
	Size gc = 0;
	for(Size i = 0; i < numberOfNodes; ++i) {
		if(SpanningForest::isClusterRep(nodes[i])) {
			if(gc < SpanningForest::clusterSize(nodes[i])) {
				gc = SpanningForest::clusterSize(nodes[i]);
			}
		}
	}
	return gc;
}

void run(
	Size numberOfNodes,
	Size initialNumberOfEdges,
	Size finalNumberOfEdges,
	Node* const nodes,
	Edge* const edges
) {
	ClusterSizeDist cdist;
	vector<Size> edgeDeletionOrder(initialNumberOfEdges);
	iota(edgeDeletionOrder.begin(), edgeDeletionOrder.end(), 0);
	shuffle(edgeDeletionOrder.begin(), edgeDeletionOrder.end(), default_random_engine(random_device()()));
	vector<Size>::iterator i = edgeDeletionOrder.begin();
	Size numberOfEdges = initialNumberOfEdges;
	Size gcSize = largestClusterSize(nodes, numberOfNodes);
	cdist[gcSize] = 1;
	while(numberOfEdges > finalNumberOfEdges) {
		disconnect(cdist, edges[*i]);
		++i;
		--numberOfEdges;
		if(cdist.rbegin()->first < gcSize) {
			gcSize = cdist.rbegin()->first;
			cout << 2.0 * numberOfEdges / numberOfNodes << "\t"
				 << 1.0 * gcSize / numberOfNodes << endl;
		}
	}
}

int main(int argc, char* argv[]) {

	if(argc == 1) {
		cout << "usage: " << argv[0] 
			 << " [number of nodes]"
			 << " [initial number of edges]"
			 << " [final number of edges]" << endl;
		return 0;
	}
	
	const Size numberOfNodes = stol(argv[1]);
	const Size initialNumberOfEdges = stod(argv[2]);
	const Size finalNumberOfEdges = stod(argv[3]);
	
	unique_ptr<Graph> g = makeERGraph(numberOfNodes, initialNumberOfEdges);
	unique_ptr<Node[]> nodes(new Node[numberOfNodes]);
	unique_ptr<Edge[]> edges(new Edge[initialNumberOfEdges]);
	initialize(*g, nodes.get(), edges.get());
	run(numberOfNodes, initialNumberOfEdges, finalNumberOfEdges, nodes.get(), edges.get());
	return 0;
}
