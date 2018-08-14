/**
 * An implementation of the decremental dynamic algorithm for mutually connected components 
 * by S. Hwang, S. Choi, D. Lee and B. Kahng (Physical Review E, 91, 022814 (2015)).
 */

#pragma once

#include "EulerTourTreeSpanningForest.hpp"
#include "DecrementalMCCModule.hpp"

namespace Snu {
namespace Cnrc {
namespace DecrementalMCC {

template<class Node, class Edge> 
using DecrementalMCCTrait = typename Snu::Cnrc::DecrementalMCCModule::DecrementalMCCModule<
	Snu::Cnrc::EulerTourTreeSpanningForest::EulerTourTreeSpanningForestNodeMixin,
	Snu::Cnrc::EulerTourTreeSpanningForest::EulerTourTreeSpanningForestEdgeMixin,
	Snu::Cnrc::EulerTourTreeSpanningForest::EulerTourTreeSpanningForestAlgorithm
>::template DecrementalMCCTrait<Node, Edge>;

template<class Node, class Edge> 
using DecrementalMCCAlgorithm = typename Snu::Cnrc::DecrementalMCCModule::DecrementalMCCModule<
	Snu::Cnrc::EulerTourTreeSpanningForest::EulerTourTreeSpanningForestNodeMixin,
	Snu::Cnrc::EulerTourTreeSpanningForest::EulerTourTreeSpanningForestEdgeMixin,
	Snu::Cnrc::EulerTourTreeSpanningForest::EulerTourTreeSpanningForestAlgorithm
>::template DecrementalMCCAlgorithm<Node, Edge>;

template<class Node, class Edge>
using DecrementalMCCNodeMixin = typename Snu::Cnrc::DecrementalMCCModule::DecrementalMCCModule<
	Snu::Cnrc::EulerTourTreeSpanningForest::EulerTourTreeSpanningForestNodeMixin,
	Snu::Cnrc::EulerTourTreeSpanningForest::EulerTourTreeSpanningForestEdgeMixin,
	Snu::Cnrc::EulerTourTreeSpanningForest::EulerTourTreeSpanningForestAlgorithm
>::template DecrementalMCCNodeMixin<Node, Edge>;

template<class Node, class Edge>
using DecrementalMCCEdgeMixin = typename Snu::Cnrc::DecrementalMCCModule::DecrementalMCCModule<
	Snu::Cnrc::EulerTourTreeSpanningForest::EulerTourTreeSpanningForestNodeMixin,
	Snu::Cnrc::EulerTourTreeSpanningForest::EulerTourTreeSpanningForestEdgeMixin,
	Snu::Cnrc::EulerTourTreeSpanningForest::EulerTourTreeSpanningForestAlgorithm
>::template DecrementalMCCEdgeMixin<Node, Edge>;

using DecrementalMCCNode = typename Snu::Cnrc::DecrementalMCCModule::DecrementalMCCModule<
	Snu::Cnrc::EulerTourTreeSpanningForest::EulerTourTreeSpanningForestNodeMixin,
	Snu::Cnrc::EulerTourTreeSpanningForest::EulerTourTreeSpanningForestEdgeMixin,
	Snu::Cnrc::EulerTourTreeSpanningForest::EulerTourTreeSpanningForestAlgorithm
>::DecrementalMCCNode;

using DecrementalMCCEdge = typename Snu::Cnrc::DecrementalMCCModule::DecrementalMCCModule<
	Snu::Cnrc::EulerTourTreeSpanningForest::EulerTourTreeSpanningForestNodeMixin,
	Snu::Cnrc::EulerTourTreeSpanningForest::EulerTourTreeSpanningForestEdgeMixin,
	Snu::Cnrc::EulerTourTreeSpanningForest::EulerTourTreeSpanningForestAlgorithm
>::DecrementalMCCEdge;

using Size = typename Snu::Cnrc::DecrementalMCCModule::DecrementalMCCModule<
	Snu::Cnrc::EulerTourTreeSpanningForest::EulerTourTreeSpanningForestNodeMixin,
	Snu::Cnrc::EulerTourTreeSpanningForest::EulerTourTreeSpanningForestEdgeMixin,
	Snu::Cnrc::EulerTourTreeSpanningForest::EulerTourTreeSpanningForestAlgorithm
>::Size;

} // end namespace DecrementalMCC
} // end namespace Cnrc
} // end namespace Snu
