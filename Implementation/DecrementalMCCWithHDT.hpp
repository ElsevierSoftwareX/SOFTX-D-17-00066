/**
 * An implementation of the decremental dynamic algorithm for mutually connected components 
 * by S. Hwang, S. Choi, D. Lee and B. Kahng (Physical Review E, 91, 022814 (2015)).
 */

#pragma once

#include "HDTSpanningForest.hpp"
#include "DecrementalMCCModule.hpp"

namespace Snu {
namespace Cnrc {
namespace DecrementalMCCWithHDT {

template<class Node, class Edge> 
using DecrementalMCCTrait = typename Snu::Cnrc::DecrementalMCCModule::DecrementalMCCModule<
	Snu::Cnrc::HDTSpanningForest::HDTSpanningForestNodeMixin,
	Snu::Cnrc::HDTSpanningForest::HDTSpanningForestEdgeMixin,
	Snu::Cnrc::HDTSpanningForest::HDTSpanningForestAlgorithm
>::template DecrementalMCCTrait<Node, Edge>;

template<class Node, class Edge> 
using DecrementalMCCAlgorithm = typename Snu::Cnrc::DecrementalMCCModule::DecrementalMCCModule<
	Snu::Cnrc::HDTSpanningForest::HDTSpanningForestNodeMixin,
	Snu::Cnrc::HDTSpanningForest::HDTSpanningForestEdgeMixin,
	Snu::Cnrc::HDTSpanningForest::HDTSpanningForestAlgorithm
>::template DecrementalMCCAlgorithm<Node, Edge>;

template<class Node, class Edge>
using DecrementalMCCNodeMixin = typename Snu::Cnrc::DecrementalMCCModule::DecrementalMCCModule<
	Snu::Cnrc::HDTSpanningForest::HDTSpanningForestNodeMixin,
	Snu::Cnrc::HDTSpanningForest::HDTSpanningForestEdgeMixin,
	Snu::Cnrc::HDTSpanningForest::HDTSpanningForestAlgorithm
>::template DecrementalMCCNodeMixin<Node, Edge>;

template<class Node, class Edge>
using DecrementalMCCEdgeMixin = typename Snu::Cnrc::DecrementalMCCModule::DecrementalMCCModule<
	Snu::Cnrc::HDTSpanningForest::HDTSpanningForestNodeMixin,
	Snu::Cnrc::HDTSpanningForest::HDTSpanningForestEdgeMixin,
	Snu::Cnrc::HDTSpanningForest::HDTSpanningForestAlgorithm
>::template DecrementalMCCEdgeMixin<Node, Edge>;

using DecrementalMCCNode = typename Snu::Cnrc::DecrementalMCCModule::DecrementalMCCModule<
	Snu::Cnrc::HDTSpanningForest::HDTSpanningForestNodeMixin,
	Snu::Cnrc::HDTSpanningForest::HDTSpanningForestEdgeMixin,
	Snu::Cnrc::HDTSpanningForest::HDTSpanningForestAlgorithm
>::DecrementalMCCNode;

using DecrementalMCCEdge = typename Snu::Cnrc::DecrementalMCCModule::DecrementalMCCModule<
	Snu::Cnrc::HDTSpanningForest::HDTSpanningForestNodeMixin,
	Snu::Cnrc::HDTSpanningForest::HDTSpanningForestEdgeMixin,
	Snu::Cnrc::HDTSpanningForest::HDTSpanningForestAlgorithm
>::DecrementalMCCEdge;

using Size = typename Snu::Cnrc::DecrementalMCCModule::DecrementalMCCModule<
	Snu::Cnrc::HDTSpanningForest::HDTSpanningForestNodeMixin,
	Snu::Cnrc::HDTSpanningForest::HDTSpanningForestEdgeMixin,
	Snu::Cnrc::HDTSpanningForest::HDTSpanningForestAlgorithm
>::Size;

} // end namespace DecrementalMCC
} // end namespace Cnrc
} // end namespace Snu
