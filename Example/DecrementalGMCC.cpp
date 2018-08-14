/**
 * This program takes three arguments [number of nodes in each network], [initial number of edges in each network]
 * and [final number of edges in each network].
 * It creates two Erdos-Renyi random networks where each has
 * [number of nodes in each network] and [initial number of edges in each network].
 * The two networks interdepend each other.
 * Each node in a network has one and only one intedependency with a node in the other network.
 * Then it prints the mean degree and the size of the largest mutually connected cluster 
 * while removing the edges one by one in a random order until [final number of edges in each network] remains in each network.
 */

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <random>

#ifdef CNRC_ENABLE_HDT
	#include "DecrementalMCCWithHDT.hpp"
	namespace MCC = Snu::Cnrc::DecrementalMCCWithHDT;
#else
	#include "DecrementalMCC.hpp"
	namespace MCC = Snu::Cnrc::DecrementalMCC;
#endif

using namespace std;

using Node = MCC::DecrementalMCCNode;
using Edge = MCC::DecrementalMCCEdge;
using NodeVector = vector<Node>;
using EdgeVector = vector<Edge>;
using MCCAlgo = MCC::DecrementalMCCAlgorithm<Node, Edge>;
using SpanningForest = MCCAlgo::SpanningForest;

bool existsEdge(const Node& n1, const Node& n2) {
	//efficient only when each node has a small degree.
	for(const Edge& e : SpanningForest::edges(n1)) {
		if(&SpanningForest::node1(e) == &n1 && &SpanningForest::node2(e) == &n2) {
			return true;
		} else if(&SpanningForest::node2(e) == &n1 && &SpanningForest::node1(e) == &n2) {
			return true;
		}
	}
	return false;
}

void makeERGraph(
	MCC::Size numberOfEdges,
	NodeVector& nodes,
	EdgeVector& edges
) {
	auto rnd = bind(
        uniform_int_distribution<unsigned int>(0, nodes.size() - 1), 
        default_random_engine(random_device()())
    );
	for(MCC::Size i = 0; i < numberOfEdges; ++i) {
		unsigned int n1 = rnd();
		unsigned int n2 = rnd();
		while(n1 == n2 || existsEdge(nodes[n1], nodes[n2])) {
			n1 = rnd();
			n2 = rnd();
		}
		SpanningForest::connect(nodes[n1], nodes[n2], edges[i]);
	}
}

void setInterdependency(NodeVector& nodes1, NodeVector& nodes2) {
	for(MCC::Size i = 0; i < nodes1.size(); ++i) {
		MCCAlgo::interconnect(nodes1[i], nodes2[i]);
	}
}

void run(
	NodeVector& nodes1,
	NodeVector& nodes2,
	EdgeVector& edges1,
	EdgeVector& edges2,
	MCC::Size finalNumberOfEdges
) {
	MCCAlgo::ClusterSizeDist cdist;
	MCCAlgo::initialize(nodes1, nodes2, cdist);
	MCC::Size gcSize = cdist.rbegin()->first;
	cout << 2.0 * edges1.size() << "\t"
		 << 1.0 * gcSize << endl;
	EdgeVector::iterator i1 = edges1.begin();
	EdgeVector::iterator i2 = edges2.begin();
	for(MCC::Size numberOfEdges = edges1.size(); 
		numberOfEdges > finalNumberOfEdges; 
		--numberOfEdges
	) {
		MCCAlgo::disconnect(*i1, cdist);
		MCCAlgo::disconnect(*i2, cdist);
		++i1;
		++i2;
		if(cdist.rbegin()->first < gcSize) {
			gcSize = cdist.rbegin()->first;
			cout << 2.0 * numberOfEdges << "\t"
				 << 1.0 * gcSize << endl;
		}
	}
}

int main(int argc, char* argv[]) {

	if(argc == 1) {
		cout << "usage: " << argv[0] 
			 << " [number of nodes in each network]"
			 << " [initial number of edges in each network]"
			 << " [final number of edges in each network]" << endl;
		return 0;
	}
	
	const MCC::Size numberOfNodes = stol(argv[1]);
	const MCC::Size initialNumberOfEdges = stol(argv[2]);
	const MCC::Size finalNumberOfEdges = stol(argv[3]);
	
	cout << setprecision(8);
	
	NodeVector nodes1(numberOfNodes);
	NodeVector nodes2(numberOfNodes);
	EdgeVector edges1(initialNumberOfEdges);
	EdgeVector edges2(initialNumberOfEdges);
	makeERGraph(initialNumberOfEdges, nodes1, edges1);
	makeERGraph(initialNumberOfEdges, nodes2, edges2);
	setInterdependency(nodes1, nodes2);
	run(nodes1, nodes2, edges1, edges2, finalNumberOfEdges);
	
	return 0;
}