#ifndef GraphAlignerBitvectorBanded_h
#define GraphAlignerBitvectorBanded_h

#include <algorithm>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <string_view>
#include "AlignmentGraph.h"
#include "NodeSlice.h"
#include "CommonUtils.h"
#include "GraphAlignerWrapper.h"
#include "ThreadReadAssertion.h"
#include "WordSlice.h"
#include "GraphAlignerBitvectorCommon.h"
#include "GraphAlignerCommon.h"
#include "ArrayPriorityQueue.h"
#include "DebugFlags.h"

template <typename LengthType, typename ScoreType, typename Word>
class GraphAlignerBitvectorBanded
{
private:
	using BV = GraphAlignerBitvectorCommon<LengthType, ScoreType, Word>;
	using Common = GraphAlignerCommon<LengthType, ScoreType, Word>;
	using AlignerGraphsizedState = typename Common::AlignerGraphsizedState;
	using Params = typename Common::Params;
	using MatrixPosition = typename Common::MatrixPosition;
	using Trace = typename Common::Trace;
	using OnewayTrace = typename Common::OnewayTrace;
	using WordSlice = typename BV::WordSlice;
	using EqVector = typename BV::EqVector;
	using EdgeWithPriority = typename Common::EdgeWithPriority;
	using DPSlice = typename BV::DPSlice;
	using DPTable = typename BV::DPTable;
	using NodeCalculationResult = typename BV::NodeCalculationResult;

public:
  /* 
   * zkn we make this public only so we can dump it to a file later. It's
   * private in all other contexts
   */
	const Params& params;

public:

	GraphAlignerBitvectorBanded(const Params& params) :
	params(params)
	{
	}


  /*
   * zkn
   * This function starts by getting a bunch of forward traces for the seed
   * hits, then it takes each forward trace and produces a backward trace from
   * the place where the forward trace started (this seems to be the equivalent
   * of a left extension, then right extension). The question is, why do we only
   * return the left extensions at the end? Don't we still want the right
   * extensions?
   * We're gonna take this one though. It appears to be the first that operates
   * on a cluster. I'm confused as to exactly how the guts of it work, but
   * evertyhing above seems to be looking more at just extracting the alignment
   * information from the traces
   */

//    // --- Add this function here ---
//     void printGraphStructure() const {
//         const AlignmentGraph& graph = params.graph;
//         std::cout << "Graph has " << graph.NodeSize() << " nodes." << std::endl;
//         for (size_t i = 0; i < graph.NodeSize(); ++i) {
//             std::cout << "Node " << i << ": ";
//             std::cout << "Name=" << graph.NodeName(i) << ", ";
//             std::cout << "Length=" << graph.NodeLength(i) << ", ";
//             std::cout << "Sequence=" << graph.NodeSequences(i, 0, graph.NodeLength(i)) << std::endl;
//             std::cout << "  InNeighbors: ";
//             for (auto n : graph.InNeighbors(i)) std::cout << n << " ";
//             std::cout << "\n  OutNeighbors: ";
//             for (auto n : graph.OutNeighbors(i)) std::cout << n << " ";
//             std::cout << std::endl;
//         }
//     }
//     // --- End function ---

void printGraphStructureForNodes(const std::vector<size_t>& nodeIds) const {
    const AlignmentGraph& graph = params.graph;
    std::ofstream out("current_slice_nodes.txt", std::ios::app); // append mode
    out << "Printing info for " << nodeIds.size() << " nodes in current slice:" << std::endl;
    for (size_t nodeId : nodeIds) {
        out << "Node " << nodeId << ": ";
        // Use the public vector for node names
        if (nodeId < graph.originalNodeName.size()) {
            out << "Name=" << graph.originalNodeName[nodeId] << ", ";
        } else {
            out << "Name=UNKNOWN, ";
        }
        out << "Length=" << graph.NodeLength(nodeId) << ", ";
        out << "Sequence=";
        for (size_t i = 0; i < graph.NodeLength(nodeId); ++i) {
            out << graph.NodeSequences(nodeId, i);
        }
        out << std::endl;
        out << "  InNeighbors: ";
        for (auto n : graph.InNeighbors(nodeId)) out << n << " ";
        out << "\n  OutNeighbors: ";
        for (auto n : graph.OutNeighbors(nodeId)) out << n << " ";
        out << std::endl;
    }
    out.close();
}


	std::vector<OnewayTrace> getMultiseedTraces(
      const std::string_view& sequence,
      const std::string_view& bwSequence,
      const std::vector<ProcessedSeedHit>& seedHits,
      AlignerGraphsizedState& reusableState,
      std::vector<ScoreType>& sliceMaxScores) const
	{

	// // Print Sequence Output
    // static size_t sequence_counter = 0;
    // sequence_counter++;
    // std::ofstream seq_out("multiseedtraces_sequence.txt", std::ios::app);
    // seq_out << "Sequence #" << sequence_counter << ": " << sequence << std::endl;
    // seq_out.close();

		std::vector<OnewayTrace> traces = getMultiseedTracesOneWay(sequence, seedHits, reusableState, sliceMaxScores, reusableState.bigraphNodeForbiddenSpans);
		std::vector<OnewayTrace> revTraces;
		std::vector<std::tuple<size_t, int, int>> reverseForbiddenNodes;
		for (auto t : reusableState.bigraphNodeForbiddenSpans)
		{
			reverseForbiddenNodes.emplace_back(int(std::get<0>(t)/2)*2 + 1 - (std::get<0>(t) % 2), sequence.size() - std::get<2>(t), sequence.size() - std::get<1>(t));
			assert(std::get<1>(reverseForbiddenNodes.back()) < std::get<2>(reverseForbiddenNodes.back()));
		}
		for (const auto& fwTrace : traces)
		{
			auto reversePos = params.graph.GetReverseDigraphPosition(fwTrace.trace[0].DPposition.node, fwTrace.trace[0].DPposition.nodeOffset);
			assert(params.graph.GetReverseDigraphPosition(reversePos.first, reversePos.second).first == fwTrace.trace[0].DPposition.node);
			assert(params.graph.GetReverseDigraphPosition(reversePos.first, reversePos.second).second == fwTrace.trace[0].DPposition.nodeOffset);

      //zkn you are flipping the read and aligning to a reversed graph, so this
      //part here is about finding the location in the flipped read from which
      //to start aligning
			assert(fwTrace.trace[0].DPposition.seqPos < sequence.size());
			size_t seqPos = sequence.size() - 1 - fwTrace.trace[0].DPposition.seqPos;
			assert(seqPos < sequence.size());

      //zkn this is a singleton vector. Seems kinda hackish
			std::vector<ProcessedSeedHit> fakeBwSeeds;
			fakeBwSeeds.emplace_back(seqPos, reversePos.first);

			std::vector<ScoreType> fakeMaxScores;
			fakeMaxScores.resize(sliceMaxScores.size(), 0);
			auto bwTraces = getMultiseedTracesOneWay(bwSequence, fakeBwSeeds, reusableState, fakeMaxScores, reverseForbiddenNodes);
			for (size_t i = 0; i < bwTraces.size(); i++)
			{
				revTraces.emplace_back(std::move(bwTraces[i]));
			}
		}
    //zkn why do we just return the reverse traces, shouldn't we return a
    //concatenation of reverse and forward traces
		return revTraces;
	}

	std::vector<OnewayTrace> getMultiseedTracesOneWay(const std::string_view& sequence, const std::vector<ProcessedSeedHit>& seedHits, AlignerGraphsizedState& reusableState, std::vector<ScoreType>& sliceMaxScores, const std::vector<std::tuple<size_t, int, int>>& forbiddenNodes) const
	{
		size_t numSlices = (sequence.size() + WordConfiguration<Word>::WordSize - 1) / WordConfiguration<Word>::WordSize;
		auto initialSlice = BV::getInitialEmptySlice();
		auto slice = getMultiseedSlices(sequence, initialSlice, numSlices, reusableState, seedHits, sliceMaxScores, forbiddenNodes);
		std::vector<OnewayTrace> results = BV::getLocalMaximaTracesFromTable(params, sequence, slice, reusableState, true, true, sliceMaxScores);
		removeDuplicateTraces(results);
		return results;
	}

	OnewayTrace getReverseTraceFromSeed(const std::string_view& sequence, int bigraphNodeId, size_t nodeOffset, int Xdropcutoff, const std::vector<std::tuple<size_t, int, int>>& forbiddenNodes, AlignerGraphsizedState& reusableState) const
	{
		size_t numSlices = (sequence.size() + WordConfiguration<Word>::WordSize - 1) / WordConfiguration<Word>::WordSize;
		auto alignmentBandwidth = BV::getInitialSliceExactPosition(params, bigraphNodeId, nodeOffset);
		auto slice = getSlices(sequence, alignmentBandwidth, numSlices, Xdropcutoff, reusableState, forbiddenNodes);
		if (slice.slices.size() <= 1)
		{
			return OnewayTrace::TraceFailed();
		}
		assert(sequence.size() <= std::numeric_limits<ScoreType>::max() - WordConfiguration<Word>::WordSize * 2);
		assert(slice.slices.back().minScore >= 0);
		assert(slice.slices.back().minScore <= (ScoreType)sequence.size() + (ScoreType)WordConfiguration<Word>::WordSize * 2);

		OnewayTrace result;
		result = BV::getReverseTraceFromTableExactEndPos(params, sequence, slice, reusableState, true, false);

		return result;
	}

	OnewayTrace getBacktraceFullStart(const std::string_view& originalSequence, int Xdropcutoff, const std::vector<std::tuple<size_t, int, int>>& forbiddenNodes, AlignerGraphsizedState& reusableState) const
	{
		assert(originalSequence.size() > 1);
		DPSlice startSlice;
		startSlice.j = -WordConfiguration<Word>::WordSize;
		startSlice.scores.addEmptyNodeMap(params.graph.NodeSize());
		startSlice.bandwidth = 1;
		startSlice.minScore = 0;
		startSlice.minScoreNode = 0;
		startSlice.minScoreNodeOffset = 0;
		startSlice.maxExactEndposScore = -params.XscoreErrorCost;
		startSlice.maxExactEndposNode = 0;
		char firstChar = originalSequence[0];
		for (size_t i = 0; i < params.graph.NodeSize(); i++)
		{
			startSlice.scores.addNodeToMap(i);
			startSlice.scores.setMinScore(i, 0);
			auto& node = startSlice.scores.node(i);
			bool match = Common::characterMatch(firstChar, params.graph.NodeSequences(i, 0));
			node.startSlice = {0, 0, match ? 0 : 1};
			node.minScore = match ? 0 : 1;
			for (size_t j = 1; j < params.graph.NodeLength(i); j++)
			{
				bool oldMatch = match;
				match = Common::characterMatch(firstChar, params.graph.NodeSequences(i, j));
				if (oldMatch && !match)
				{
					node.HP[j / params.graph.SPLIT_NODE_SIZE] |= ((Word)1) << (j % params.graph.SPLIT_NODE_SIZE);
				}
				else if (match && !oldMatch)
				{
					node.HN[j / params.graph.SPLIT_NODE_SIZE] |= ((Word)1) << (j % params.graph.SPLIT_NODE_SIZE);
				}
				if (match) node.minScore = 0;
			}
			if (node.minScore == 0)
			{
				startSlice.maxExactEndposScore = 0;
				startSlice.maxExactEndposNode = i;
			}
			node.endSlice = {0, 0, match ? 0 : 1};
			node.exists = true;
		}
		std::string_view alignableSequence { originalSequence.data()+1, originalSequence.size() - 1 };
		assert(alignableSequence.size() > 0);
		size_t numSlices = (alignableSequence.size() + WordConfiguration<Word>::WordSize - 1) / WordConfiguration<Word>::WordSize;
		auto slice = getSlices(alignableSequence, startSlice, numSlices, Xdropcutoff, reusableState, forbiddenNodes);
		if (slice.slices.size() <= 1)
		{
			return OnewayTrace::TraceFailed();
		}

		OnewayTrace result;
		result = BV::getReverseTraceFromTableExactEndPos(params, alignableSequence, slice, reusableState, true, false);
		for (size_t i = 0; i < result.trace.size(); i++)
		{
			result.trace[i].DPposition.seqPos += 1;
		}
		std::reverse(result.trace.begin(), result.trace.end());
		result.trace[0].sequenceCharacter = originalSequence[0];
		assert(result.trace[0].DPposition.seqPos == 0);
		return result;
	}

private:

	int traceAlignmentScore(const OnewayTrace& trace) const
	{
		assert(trace.trace[0].DPposition.seqPos >= trace.trace.back().DPposition.seqPos);
		ScoreType score = (trace.trace[0].DPposition.seqPos - trace.trace.back().DPposition.seqPos + 1)*100;
		score -= (ScoreType)trace.score * params.XscoreErrorCost;
		return score;
	}

	void removeDifferentTracePrefix(OnewayTrace& removeFrom, const OnewayTrace& removeWith) const
	{
		size_t lastEqualOffset = 0;
		assert(removeFrom.trace.size() > 0);
		assert(removeWith.trace.size() > 0);
		assert(removeFrom.trace.back() == removeWith.trace.back());
		while (lastEqualOffset+1 < removeFrom.trace.size() && lastEqualOffset+1 < removeWith.trace.size() && removeFrom.trace[removeFrom.trace.size()-1-(lastEqualOffset+1)] == removeWith.trace[removeWith.trace.size()-1-(lastEqualOffset+1)])
		{
			lastEqualOffset += 1;
		}
		assert(lastEqualOffset < removeFrom.trace.size());
		assert(lastEqualOffset > 0);
		removeFrom.trace.erase(removeFrom.trace.begin(), removeFrom.trace.end() - lastEqualOffset);
	}

	void removeDuplicateTraces(std::vector<OnewayTrace>& traces) const
	{
		if (traces.size() == 0) return;
		std::sort(traces.begin(), traces.end(), [this](const OnewayTrace& left, const OnewayTrace& right) {
			return traceAlignmentScore(left) > traceAlignmentScore(right);
		});
		for (size_t i = traces.size()-1; i > 0; i--)
		{
			for (size_t j = 0; j < i; j++)
			{
				if (traces[i].trace.back().DPposition == traces[j].trace.back().DPposition)
				{
					// times 100 because int scores are 1/100ths
					if (params.clipAmbiguousEnds >= 0 && traceAlignmentScore(traces[i]) >= traceAlignmentScore(traces[j]) - params.clipAmbiguousEnds * 100)
					{
						removeDifferentTracePrefix(traces[j], traces[i]);
					}
					std::swap(traces[i], traces.back());
					traces.pop_back();
					break;
				}
			}
		}
	}

#ifdef EXTRACORRECTNESSASSERTIONS
	template <bool HasVectorMap, bool PreviousHasVectorMap>
	void checkNodeBoundaryCorrectness(const NodeSlice<LengthType, ScoreType, Word, HasVectorMap>& currentSlice, const NodeSlice<LengthType, ScoreType, Word, PreviousHasVectorMap>& previousSlice, const std::string_view& sequence, size_t j, ScoreType maxScore, ScoreType previousMaxScore, const std::vector<bool>& hasSeedStart, const WordSlice seedstartSlice, const WordSlice fakeSlice) const
	{
		assert(previousMaxScore <= maxScore || seedstartSlice.getScoreBeforeStart() <= maxScore);
		for (auto pair : currentSlice)
		{
			auto node = pair.first;
			WordSlice extraSlice = hasSeedStart[node] ? seedstartSlice : fakeSlice;
			if (previousSlice.hasNode(node) && previousSlice.node(node).exists)
			{
				assert(pair.second.startSlice.getScoreBeforeStart() <= previousSlice.node(node).startSlice.scoreEnd);
			}
			bool eq = Common::characterMatch(sequence[j], params.graph.NodeSequences(node, 0));
			ScoreType foundMinScore = std::numeric_limits<ScoreType>::max();
			foundMinScore = std::min(foundMinScore, currentSlice.node(node).startSlice.getScoreBeforeStart()+1);
			if (previousSlice.hasNode(node) && previousSlice.node(node).exists && previousSlice.node(node).minScore <= previousMaxScore)
			{
				foundMinScore = std::min(foundMinScore, previousSlice.node(node).startSlice.scoreEnd+1);
			}
			for (auto neighbor : params.graph.InNeighbors(node))
			{
				if (currentSlice.hasNode(neighbor) && currentSlice.node(neighbor).exists)
				{
					foundMinScore = std::min(foundMinScore, currentSlice.node(neighbor).endSlice.getValue(0)+1);
					foundMinScore = std::min(foundMinScore, currentSlice.node(neighbor).endSlice.getScoreBeforeStart() + (eq?0:1));
				}
				if (previousSlice.hasNode(neighbor) && previousSlice.node(neighbor).exists && previousSlice.node(neighbor).minScore <= previousMaxScore)
				{
					foundMinScore = std::min(foundMinScore, previousSlice.node(neighbor).endSlice.scoreEnd + (eq ? 0 : 1));
				}
			}
			if (pair.second.startSlice.getValue(0) <= maxScore || foundMinScore <= maxScore)
			{
				assert(foundMinScore != std::numeric_limits<ScoreType>::max());
				// todo figure out why breaks
				assert(pair.second.startSlice.getValue(0) == foundMinScore || pair.second.startSlice.getValue(0) == extraSlice.getValue(0));
			}
			bool hasMatch = Common::characterMatch(sequence[0], params.graph.NodeSequences(node, 0));
			for (size_t i = 1; i < WordConfiguration<Word>::WordSize; i++)
			{
				if (j+i >= sequence.size()) break;
				eq = Common::characterMatch(sequence[j+i], params.graph.NodeSequences(node, 0));
				if (eq) hasMatch = true;
				ScoreType foundMinScore = pair.second.startSlice.getValue(i-1)+1;
				for (auto neighbor : params.graph.InNeighbors(node))
				{
					if (!currentSlice.hasNode(neighbor)) continue;
					if (!currentSlice.node(neighbor).exists) continue;
					foundMinScore = std::min(foundMinScore, currentSlice.node(neighbor).endSlice.getValue(i)+1);
					foundMinScore = std::min(foundMinScore, currentSlice.node(neighbor).endSlice.getValue(i-1) + (eq ? 0 : 1));
				}
				if (pair.second.startSlice.getValue(i) <= maxScore || foundMinScore <= maxScore)
				{
					// todo figure out why breaks
					assert(pair.second.startSlice.getValue(i) == foundMinScore || pair.second.startSlice.getValue(i) == extraSlice.getValue(i));
				}
			}
		}
	}
#endif

	template <bool HasVectorMap, bool PreviousHasVectorMap, typename PriorityQueue>
	void addSeedHitToScoresAndQueue(const ProcessedSeedHit& seedHit, NodeSlice<LengthType, ScoreType, Word, HasVectorMap>& currentSlice, const NodeSlice<LengthType, ScoreType, Word, PreviousHasVectorMap>& previousSlice, std::vector<bool>& currentBand, const std::vector<bool>& previousBand, PriorityQueue& calculableQueue, const WordSlice extraSlice, phmap::flat_hash_map<size_t, ScoreType>& nodeMaxExactEndposScore, bool storeNodeExactEndposScores) const
	{
#ifdef SLICEVERBOSE
		std::cerr << " " << seedHit.alignmentGraphNodeId << "(" << params.graph.BigraphNodeID(seedHit.alignmentGraphNodeId) << ")";
#endif
		size_t node = seedHit.alignmentGraphNodeId;
		assert(calculableQueue.IsComponentPriorityQueue());
		if (calculableQueue.IsComponentPriorityQueue())
		{
			EdgeWithPriority insertEdge { node, extraSlice.getValue(0), extraSlice, true };
			insertEdge.forceCalculation = true;
			calculableQueue.insert(params.graph.ComponentNumber(node), extraSlice.getValue(0), insertEdge);
		}
	}

	template <bool HasVectorMap, bool PreviousHasVectorMap, typename PriorityQueue>
	__attribute__((noinline))
	NodeCalculationResult calculateSlice(const std::string_view& sequence, const size_t j, NodeSlice<LengthType, ScoreType, Word, HasVectorMap>& currentSlice, 
		const NodeSlice<LengthType, ScoreType, Word, PreviousHasVectorMap>& previousSlice, std::vector<bool>& currentBand, const std::vector<bool>& previousBand, 
		PriorityQueue& calculableQueue, ScoreType previousQuitScore, ScoreType bandwidth, ScoreType previousMinScore, const std::vector<ProcessedSeedHit>& seedHits, 
		size_t seedhitStart, size_t seedhitEnd, const WordSlice seedstartSlice, std::vector<bool>& hasSeedStart, std::unordered_set<size_t>& seedstartNodes, 
		phmap::flat_hash_map<size_t, ScoreType>& nodeMaxExactEndposScore, bool storeNodeExactEndposScores, const std::vector<bool>& allowedBigraphNodesThisSlice) const
	{



		// std::ofstream calc_out("calculateslice_sequence.txt", std::ios::app);
		// calc_out << "Offset: " << j << "Sequence: " << sequence << std::endl;
		// calc_out.close();

		// --- Debug logging block ---
		if (debugTop == true) {
		if (enableCalculateSliceDebug) {
			calculateSliceIteration++;
			std::ofstream dbg("GbvCallTrace.log", std::ios::app);
			dbg << "calculateSlice call #" << calculateSliceIteration
				<< " | j=" << j
				<< " | calculableQueue.size()=" << calculableQueue.size()
				<< " | sequence[0:64]=";
			for (size_t idx = 0; idx < std::min<size_t>(64, sequence.size()); ++idx) {
				dbg << sequence[idx];
			}
			dbg << std::endl;
			dbg.close();
		}

		if (
			calculateSliceIteration >= 3 &&
			getNextSliceIteration >= 3 &&
			calculateNodeClipPreciseIteration >= 3 &&
			calculateNodeInnerIteration >= 3 &&
			getScoreBeforeStartIteration >= 3 &&
			mergeTwoSlices2InputIteration >= 3 &&
			mergeTwoSlices4InputIteration >= 3 &&
			differenceMasksBitTwiddleIteration >= 3 &&
			flattenWordSliceIteration >= 3
		) {
			enableCalculateSliceDebug = false;
			enableGetNextSliceDebug = false;
			enableCalculateNodeClipPreciseDebug = false;
			enableCalculateNodeInnerDebug = false;
			enableGetScoreBeforeStartDebug = false;
			enableMergeTwoSlices2InputDebug = false;
			enableMergeTwoSlices4InputDebug = false;
			enableDifferenceMasksBitTwiddleDebug = false;
			enableFlattenWordSliceDebug = false;
		}
	}
		// --- End debug logging block ---

		

		if (previousMinScore == std::numeric_limits<ScoreType>::max() - bandwidth - 1)
		{
			assert(seedstartSlice.scoreEnd != std::numeric_limits<ScoreType>::max());
			previousMinScore = seedstartSlice.getScoreBeforeStart();
			previousQuitScore = seedstartSlice.getScoreBeforeStart() + bandwidth;
		}
		if (seedstartSlice.scoreEnd != std::numeric_limits<ScoreType>::max() && previousMinScore > seedstartSlice.scoreEnd)
		{
			previousMinScore = seedstartSlice.scoreEnd;
			previousQuitScore = std::max(previousQuitScore, seedstartSlice.scoreEnd + bandwidth);
		}
		if (j > 0)
		{
			bool hasPreviousNode = false;
			for (auto node : previousSlice)
			{
				if (!(allowedBigraphNodesThisSlice[params.graph.BigraphNodeID(node.first)])) continue;
				hasPreviousNode = true;
				if (node.second.minScore > previousQuitScore) continue;
				if (params.graph.Linearizable(node.first))
				{
					auto neighbor = params.graph.InNeighbors(node.first)[0];
					if (previousBand[neighbor] && previousSlice.node(neighbor).endSlice.scoreEnd < previousQuitScore && previousSlice.node(neighbor).minScore < previousQuitScore)
					{
						//linear area, no need to add the later node into the queue
						//because calculating the earlier node will guarantee that the later node will get added
						continue;
					}
				}
				hasPreviousNode = true;
			}
			if (!hasPreviousNode)
			{
				previousQuitScore = std::max(previousQuitScore, seedstartSlice.scoreEnd + bandwidth);
			}
		}
		double averageErrorRate = 0;
		if (j > 0)
		{
			averageErrorRate = (double)previousMinScore / (double)j;
		}
		double priorityMismatchPenalty = WordConfiguration<Word>::WordSize;
		if (averageErrorRate > 1.0/(double)WordConfiguration<Word>::WordSize)
		{
			priorityMismatchPenalty = 1.0 / averageErrorRate;
		}
		assert(priorityMismatchPenalty >= 0);
		assert(priorityMismatchPenalty <= WordConfiguration<Word>::WordSize);
		NodeCalculationResult result;
		result.minScore = std::numeric_limits<ScoreType>::max() - bandwidth - 1;
		result.minScoreNode = std::numeric_limits<LengthType>::max();
		result.minScoreNodeOffset = std::numeric_limits<LengthType>::max();
		result.maxExactEndposNode = std::numeric_limits<LengthType>::min();
		result.maxExactEndposScore = std::numeric_limits<ScoreType>::min();
		result.cellsProcessed = 0;
#ifdef SLICEVERBOSE
		result.nodesProcessed = 0;
#endif

		EqVector EqV = BV::getEqVector(sequence, j);

		assert(previousSlice.size() > 0 || seedhitStart != std::numeric_limits<size_t>::max());
		ScoreType zeroScore = previousMinScore*priorityMismatchPenalty - j - 64;
		if (j == 0)
		{
			for (auto node : previousSlice)
			{
				if (!(allowedBigraphNodesThisSlice[params.graph.BigraphNodeID(node.first)])) continue;
				assert(node.second.minScore <= previousQuitScore);
				WordSlice startSlice = BV::getSourceSliceFromScore(node.second.startSlice.scoreEnd);
				if (calculableQueue.IsComponentPriorityQueue())
				{
					calculableQueue.insert(params.graph.ComponentNumber(node.first), node.second.minScore, EdgeWithPriority { node.first, node.second.minScore - previousMinScore, startSlice, true });
				}
				else
				{
					calculableQueue.insert(node.second.minScore*priorityMismatchPenalty - j - zeroScore, EdgeWithPriority { node.first, node.second.minScore - previousMinScore, startSlice, true });
				}
			}
		}
		else
		{
			for (auto node : previousSlice)
			{
				if (!(allowedBigraphNodesThisSlice[params.graph.BigraphNodeID(node.first)])) continue;
				assert(node.second.exists);
				assert(node.second.minScore <= node.second.startSlice.scoreEnd);
				assert(node.second.minScore <= node.second.endSlice.scoreEnd);
				if (node.second.minScore > previousQuitScore) continue;
				if (params.graph.Linearizable(node.first))
				{
					auto neighbor = params.graph.InNeighbors(node.first)[0];
				 	if (previousBand[neighbor] && previousSlice.node(neighbor).endSlice.scoreEnd < previousQuitScore && previousSlice.node(neighbor).minScore < previousQuitScore)
				 	{
				 		//linear area, no need to add the later node into the queue 
				 		//because calculating the earlier node will guarantee that the later node will get added
				 		continue;
				 	}
				}
				WordSlice startSlice = BV::getSourceSliceFromScore(node.second.startSlice.scoreEnd);
				if (calculableQueue.IsComponentPriorityQueue())
				{
					calculableQueue.insert(params.graph.ComponentNumber(node.first), node.second.minScore, EdgeWithPriority { node.first, node.second.minScore - previousMinScore, startSlice, true });
				}
				else
				{
					calculableQueue.insert(node.second.minScore*priorityMismatchPenalty - j - zeroScore, EdgeWithPriority { node.first, node.second.minScore - previousMinScore, startSlice, true });
				}
			}
		}
		assert(calculableQueue.size() > 0 || seedhitStart != std::numeric_limits<size_t>::max());

#ifdef SLICEVERBOSE
		std::cerr << "seeds";
#endif

		std::vector<size_t> clearSeedStarts;
		if (seedstartSlice.getScoreBeforeStart() < previousQuitScore + WordConfiguration<Word>::WordSize)
		{
#ifdef SLICEVERBOSE
			std::cerr << " " << seedhitEnd - seedhitStart << ":";
#endif
			for (size_t i = seedhitStart; i < seedhitEnd; i++)
			{
				if (hasSeedStart[seedHits[i].alignmentGraphNodeId]) continue;
				assert(allowedBigraphNodesThisSlice[params.graph.BigraphNodeID(seedHits[i].alignmentGraphNodeId)]);
				addSeedHitToScoresAndQueue(seedHits[i], currentSlice, previousSlice, currentBand, previousBand, calculableQueue, seedstartSlice, nodeMaxExactEndposScore, storeNodeExactEndposScores);
				assert(!hasSeedStart[seedHits[i].alignmentGraphNodeId]);
				hasSeedStart[seedHits[i].alignmentGraphNodeId] = true;
				clearSeedStarts.push_back(seedHits[i].alignmentGraphNodeId);
				seedstartNodes.insert(seedHits[i].alignmentGraphNodeId);
			}
#ifdef SLICEVERBOSE
			std::cerr << std::endl;
#endif
		}
		else
		{
#ifdef SLICEVERBOSE
			std::cerr << " not added because of score" << std::endl;
#endif
		}

		WordSlice fakeSlice { WordConfiguration<Word>::AllZeros, WordConfiguration<Word>::AllZeros, std::numeric_limits<ScoreType>::max() };
		ScoreType currentMinScoreAtEndRow = result.minScore;

		// Main Queue Here 

		while (calculableQueue.size() > 0)
		{

		// --- Debug logging block ---
		if (debugTop == true) {
		if (enableCalculateSliceDebug) {
			std::ofstream dbg("GbvCallTrace.log", std::ios::app);
			dbg << "calculateSlice Queue Start call #" << calculateSliceIteration
			<< " | calculableQueue.size()=" << calculableQueue.size()
			<< std::endl;
			dbg.close();
			}

		if (
			calculateSliceIteration >= 3 &&
			getNextSliceIteration >= 3 &&
			calculateNodeClipPreciseIteration >= 3 &&
			calculateNodeInnerIteration >= 3 &&
			getScoreBeforeStartIteration >= 3 &&
			mergeTwoSlices2InputIteration >= 3 &&
			mergeTwoSlices4InputIteration >= 3 &&
			differenceMasksBitTwiddleIteration >= 3 &&
			flattenWordSliceIteration >= 3
		) {
			enableCalculateSliceDebug = false;
			enableGetNextSliceDebug = false;
			enableCalculateNodeClipPreciseDebug = false;
			enableCalculateNodeInnerDebug = false;
			enableGetScoreBeforeStartDebug = false;
			enableMergeTwoSlices2InputDebug = false;
			enableMergeTwoSlices4InputDebug = false;
			enableDifferenceMasksBitTwiddleDebug = false;
			enableFlattenWordSliceDebug = false;
		}
	}
		// --- End debug logging block ---

			auto pair = calculableQueue.top();
			if (!calculableQueue.IsComponentPriorityQueue())
			{
				if (pair.priority > currentMinScoreAtEndRow + bandwidth) break;
			}
			if (calculableQueue.extraSize(pair.target) == 0)
			{
				calculableQueue.pop();
				continue;
			}
			auto i = pair.target;
			//NEW CODE ADDED FOR GRAPH DEBUG
			printGraphStructureForNodes(std::vector<size_t>{i});
			bool firstCalc = false;
			assert(allowedBigraphNodesThisSlice[params.graph.BigraphNodeID(i)]);
			if (!currentBand[i])
			{
				assert(!currentSlice.hasNode(i));
				currentSlice.addNode(i);
				currentBand[i] = true;
				firstCalc = true;
				if (hasSeedStart[i])
				{
					currentSlice.addNodeToMap(i);
					currentSlice.setMinScore(i, seedstartSlice.scoreEnd);
					auto& nodeScores = currentSlice.node(i);
					nodeScores.startSlice = seedstartSlice;
					nodeScores.endSlice = seedstartSlice;
					nodeScores.exists = true;
				}
			}
			assert(currentBand[i]);
			const std::vector<EdgeWithPriority>* extras;
			extras = &calculableQueue.getExtras(i);
			auto& thisNode = currentSlice.node(i);
			auto oldEnd = thisNode.endSlice;
			if (!thisNode.exists) oldEnd = { 0, 0, std::numeric_limits<ScoreType>::max() };
			WordSlice extraSlice = hasSeedStart[i] ? seedstartSlice : fakeSlice;
			if (extraSlice.scoreEnd != std::numeric_limits<ScoreType>::max()) oldEnd = oldEnd.mergeWith(extraSlice);
			typename NodeSlice<LengthType, ScoreType, Word, HasVectorMap>::NodeSliceMapItem previousThisNode;

			if (previousBand[i])
			{
				previousThisNode = previousSlice.node(i);
				assert(previousThisNode.exists);
			}
			else
			{
				for (size_t chunk = 0; chunk < previousThisNode.NUM_CHUNKS; chunk++)
				{
					previousThisNode.HP[chunk] = WordConfiguration<Word>::AllOnes;
					previousThisNode.HN[chunk] = WordConfiguration<Word>::AllZeros;
				}
				previousThisNode.exists = false;
			}
			NodeCalculationResult nodeCalc;
#ifdef SLICEVERBOSE
			std::cerr << "calculate node " << i << " (" << params.graph.BigraphNodeID(i) << ")" << std::endl;
#endif
			if (i < params.graph.FirstAmbiguous())
			{
				nodeCalc = BV::calculateNodeClipPrecise(params, i, thisNode, EqV, previousThisNode, *extras, previousBand, params.graph.NodeChunks(i), extraSlice, j);
				assert(nodeCalc.maxExactEndposScore != std::numeric_limits<ScoreType>::min());
			}
			else
			{
				nodeCalc = BV::calculateNodeClipPrecise(params, i, thisNode, EqV, previousThisNode, *extras, previousBand, params.graph.AmbiguousNodeChunks(i), extraSlice, j);
				assert(nodeCalc.maxExactEndposScore != std::numeric_limits<ScoreType>::min());
			}
			assert(nodeCalc.minScore != std::numeric_limits<ScoreType>::max());
			calculableQueue.pop();
			if (!calculableQueue.IsComponentPriorityQueue())
			{
				calculableQueue.removeExtras(i);
			}
			if (storeNodeExactEndposScores)
			{
				nodeMaxExactEndposScore[i] = std::max(nodeCalc.maxExactEndposScore, nodeMaxExactEndposScore[i]);
			}
			// todo fix
			// assert(nodeCalc.minScore <= previousQuitScore + bandwidth + params.graph.SPLIT_NODE_SIZE + WordConfiguration<Word>::WordSize);
			assert(nodeCalc.minScore <= (ScoreType)j + bandwidth + (ScoreType)WordConfiguration<Word>::WordSize + (ScoreType)WordConfiguration<Word>::WordSize);
			currentMinScoreAtEndRow = std::min(currentMinScoreAtEndRow, nodeCalc.minScore);
			currentSlice.setMinScoreIfSmaller(i, nodeCalc.minScore);
#ifdef SLICEVERBOSE
			volatile size_t firstslices = currentSlice.node(i).firstSlicesCalcedWhenCalced;
			volatile size_t calcedslices = currentSlice.node(i).slicesCalcedWhenCalced;
			if (currentSlice.node(i).firstSlicesCalcedWhenCalced == std::numeric_limits<size_t>::max()) currentSlice.node(i).firstSlicesCalcedWhenCalced = result.cellsProcessed;
			if (currentSlice.node(i).slicesCalcedWhenCalced != std::numeric_limits<size_t>::max()) assert(currentSlice.node(i).slicesCalcedWhenCalced < result.cellsProcessed);
			currentSlice.node(i).slicesCalcedWhenCalced = result.cellsProcessed;
			assert(currentSlice.node(i).firstSlicesCalcedWhenCalced <= currentSlice.node(i).slicesCalcedWhenCalced);
#endif
			auto newEnd = thisNode.endSlice;

			if (firstCalc || newEnd.scoreEnd != oldEnd.scoreEnd || newEnd.VP != oldEnd.VP || newEnd.VN != oldEnd.VN)
			{
				ScoreType newEndMinScore = newEnd.changedMinScore(oldEnd);
				if (firstCalc) newEndMinScore = newEnd.getMinScore();
				// assert(newEndMinScore >= previousMinScore || newEndMinScore >= seedstartSlice.getScoreBeforeStart());
				assert(newEndMinScore != std::numeric_limits<ScoreType>::max());
				if (newEndMinScore <= currentMinScoreAtEndRow + bandwidth)
				{
					for (auto neighbor : params.graph.OutNeighbors(i))
					{
						if (!allowedBigraphNodesThisSlice[params.graph.BigraphNodeID(neighbor)]) continue;
						if (calculableQueue.IsComponentPriorityQueue())
						{
							calculableQueue.insert(params.graph.ComponentNumber(neighbor), newEndMinScore, EdgeWithPriority { neighbor, newEndMinScore - previousMinScore, newEnd, false });
						}
						else
						{
							ScoreType newEndPriorityScore = newEnd.getChangedPriorityScore(oldEnd, j, priorityMismatchPenalty);
							assert(newEndPriorityScore != std::numeric_limits<ScoreType>::max());
							assert(newEndPriorityScore >= zeroScore);
							calculableQueue.insert(newEndPriorityScore - zeroScore, EdgeWithPriority { neighbor, newEndMinScore - previousMinScore, newEnd, false });
						}
					}
				}
			}
			if (nodeCalc.minScore < result.minScore)
			{
				result.minScore = nodeCalc.minScore;
				result.minScoreNode = nodeCalc.minScoreNode;
				result.minScoreNodeOffset = nodeCalc.minScoreNodeOffset;
			}
			if (nodeCalc.maxExactEndposScore > result.maxExactEndposScore)
			{
				result.maxExactEndposScore = nodeCalc.maxExactEndposScore;
				result.maxExactEndposNode = nodeCalc.maxExactEndposNode;
			}
			assert(currentSlice.node(i).minScore <= currentSlice.node(i).startSlice.scoreEnd);
			assert(currentSlice.node(i).minScore <= currentSlice.node(i).endSlice.scoreEnd);
			assert(result.minScore == currentMinScoreAtEndRow);
			result.cellsProcessed += nodeCalc.cellsProcessed;
			assert(nodeCalc.cellsProcessed > 0);
#ifdef SLICEVERBOSE
			result.nodesProcessed++;
#endif
			if (result.cellsProcessed > params.maxCellsPerSlice) break;

		// --- Debug logging block ---
		if (debugTop == true) {
		if (enableCalculateSliceDebug) {
			std::ofstream dbg("GbvCallTrace.log", std::ios::app);
			dbg << "calculateSlice Queue End call #" << calculateSliceIteration
			<< " | calculableQueue.size()=" << calculableQueue.size()
			<< std::endl;
			dbg.close();
			}

		if (
			calculateSliceIteration >= 3 &&
			getNextSliceIteration >= 3 &&
			calculateNodeClipPreciseIteration >= 3 &&
			calculateNodeInnerIteration >= 3 &&
			getScoreBeforeStartIteration >= 3 &&
			mergeTwoSlices2InputIteration >= 3 &&
			mergeTwoSlices4InputIteration >= 3 &&
			differenceMasksBitTwiddleIteration >= 3 &&
			flattenWordSliceIteration >= 3
		) {
			enableCalculateSliceDebug = false;
			enableGetNextSliceDebug = false;
			enableCalculateNodeClipPreciseDebug = false;
			enableCalculateNodeInnerDebug = false;
			enableGetScoreBeforeStartDebug = false;
			enableMergeTwoSlices2InputDebug = false;
			enableMergeTwoSlices4InputDebug = false;
			enableDifferenceMasksBitTwiddleDebug = false;
			enableFlattenWordSliceDebug = false;
		}
	}
		// --- End debug logging block ---

		} //end while

		// END OF QUEUE

#ifdef EXTRACORRECTNESSASSERTIONS
		checkNodeBoundaryCorrectness<HasVectorMap, PreviousHasVectorMap>(currentSlice, previousSlice, sequence, j, currentMinScoreAtEndRow + bandwidth, previousQuitScore, hasSeedStart, seedstartSlice, fakeSlice);
#endif
		for (auto node : currentSlice)
		{
			assert(node.second.exists);
			assert(node.second.minScore <= node.second.startSlice.scoreEnd);
			assert(node.second.minScore <= node.second.endSlice.scoreEnd);
		}

		for (auto node : clearSeedStarts)
		{
			assert(hasSeedStart[node]);
			hasSeedStart[node] = false;
		}

		assert(result.minScoreNode != std::numeric_limits<LengthType>::max() || (seedhitStart == seedhitEnd && seedhitStart != std::numeric_limits<size_t>::max()));

#ifdef SLICEVERBOSE
		std::cerr << "prefilternodes " << currentSlice.size() << " ";
#endif

		calculableQueue.clear();

		return result;
	}

	template <typename PriorityQueue>
	void fillDPSlice(const std::string_view& sequence, DPSlice& slice, const DPSlice& previousSlice, const std::vector<bool>& previousBand, std::vector<bool>& currentBand, PriorityQueue& calculableQueue, ScoreType bandwidth, const std::vector<ProcessedSeedHit>& seedHits, size_t seedhitStart, size_t seedhitEnd, const WordSlice extraSlice, std::vector<bool>& hasSeedStart, bool storeNodeExactEndposScores, const std::vector<bool>& allowedBigraphNodesThisSlice) const
	{
		NodeCalculationResult sliceResult;
		assert((ScoreType)previousSlice.bandwidth < std::numeric_limits<ScoreType>::max());
		assert((ScoreType)previousSlice.bandwidth >= 0);
		assert(previousSlice.minScore < std::numeric_limits<ScoreType>::max() - (ScoreType)previousSlice.bandwidth);
		if (slice.scoresVectorMap.hasVectorMapCurrently())
		{
			if (previousSlice.scoresVectorMap.hasVectorMapCurrently())
			{
				sliceResult = calculateSlice<true, true>(sequence, slice.j, slice.scoresVectorMap, previousSlice.scoresVectorMap, currentBand, previousBand, calculableQueue, previousSlice.minScore + previousSlice.bandwidth, bandwidth, previousSlice.minScore, seedHits, seedhitStart, seedhitEnd, extraSlice, hasSeedStart, slice.seedstartNodes, slice.nodeMaxExactEndposScore, storeNodeExactEndposScores, allowedBigraphNodesThisSlice);
			}
			else
			{
				sliceResult = calculateSlice<true, false>(sequence, slice.j, slice.scoresVectorMap, previousSlice.scores, currentBand, previousBand, calculableQueue, previousSlice.minScore + previousSlice.bandwidth, bandwidth, previousSlice.minScore, seedHits, seedhitStart, seedhitEnd, extraSlice, hasSeedStart, slice.seedstartNodes, slice.nodeMaxExactEndposScore, storeNodeExactEndposScores, allowedBigraphNodesThisSlice);
			}
			slice.scores = slice.scoresVectorMap.getMapSlice();
		}
		else
		{
			assert(!previousSlice.scoresVectorMap.hasVectorMapCurrently());
			sliceResult = calculateSlice<false, false>(sequence, slice.j, slice.scores, previousSlice.scores, currentBand, previousBand, calculableQueue, previousSlice.minScore + previousSlice.bandwidth, bandwidth, previousSlice.minScore, seedHits, seedhitStart, seedhitEnd, extraSlice, hasSeedStart, slice.seedstartNodes, slice.nodeMaxExactEndposScore, storeNodeExactEndposScores, allowedBigraphNodesThisSlice);
		}
		slice.cellsProcessed = sliceResult.cellsProcessed;
		slice.minScoreNode = sliceResult.minScoreNode;
		slice.minScoreNodeOffset = sliceResult.minScoreNodeOffset;
		slice.minScore = sliceResult.minScore;
		slice.maxExactEndposScore = sliceResult.maxExactEndposScore;
		slice.maxExactEndposNode = sliceResult.maxExactEndposNode;
		assert((seedHits.size() != 0 && seedhitEnd == seedhitStart) || sliceResult.minScore <= (ScoreType)slice.j + (ScoreType)WordConfiguration<Word>::WordSize);
		assert((seedHits.size() != 0 && seedhitEnd == seedhitStart) || sliceResult.maxExactEndposScore != std::numeric_limits<ScoreType>::min());
		assert((seedHits.size() != 0 && seedhitEnd == seedhitStart) || (double)sliceResult.maxExactEndposScore >= -((double)slice.j + WordConfiguration<Word>::WordSize) * params.XscoreErrorCost);
		assert((seedHits.size() != 0 && seedhitEnd == seedhitStart) || (double)sliceResult.maxExactEndposScore <= ((double)slice.j + WordConfiguration<Word>::WordSize) * params.XscoreErrorCost);
		assert((seedHits.size() != 0 && seedhitEnd == seedhitStart) || (double)slice.maxExactEndposScore <= ((double)slice.j + WordConfiguration<Word>::WordSize) * params.XscoreErrorCost);
		assert((seedHits.size() != 0 && seedhitEnd == seedhitStart) || (double)slice.maxExactEndposScore >= -((double)slice.j + WordConfiguration<Word>::WordSize) * params.XscoreErrorCost);
		assert(slice.minScore >= previousSlice.minScore || slice.minScore >= extraSlice.getValue(0) || previousSlice.minScore == std::numeric_limits<ScoreType>::max() - bandwidth - 1);
		slice.bandwidth = bandwidth;
#ifdef SLICEVERBOSE
		slice.nodesProcessed = sliceResult.nodesProcessed;
		for (auto node : slice.scores)
		{
			if (currentBand[node.first])
			{
				slice.numCells += params.graph.NodeLength(node.first);
			}
		}
#endif
	}

	template <typename PriorityQueue>
	DPSlice pickMethodAndExtendFill(const std::string_view& sequence, const DPSlice& previous, const std::vector<bool>& previousBand, std::vector<bool>& currentBand, PriorityQueue& calculableQueue, ScoreType bandwidth, const std::vector<ProcessedSeedHit>& seedHits, size_t seedhitStart, size_t seedhitEnd, const WordSlice extraSlice, std::vector<bool>& hasSeedStart, bool storeNodeExactEndposScores, const std::vector<bool>& allowedBigraphNodesThisSlice) const
	{
		DPSlice bandTest;
		bandTest.scores.addEmptyNodeMap(previous.scores.size());
		bandTest.j = previous.j + WordConfiguration<Word>::WordSize;
		fillDPSlice(sequence, bandTest, previous, previousBand, currentBand, calculableQueue, bandwidth, seedHits, seedhitStart, seedhitEnd, extraSlice, hasSeedStart, storeNodeExactEndposScores, allowedBigraphNodesThisSlice);
		return bandTest;
	}

	DPTable getSlices(const std::string_view& sequence, const DPSlice& initialSlice, size_t numSlices, int Xdropcutoff, AlignerGraphsizedState& reusableState, const std::vector<std::tuple<size_t, int, int>>& forbiddenNodes) const
	{
		return getXdropSlices(sequence, initialSlice, numSlices, Xdropcutoff, forbiddenNodes, reusableState);
	}

	DPTable getXdropSlices(const std::string_view& sequence, const DPSlice& initialSlice, size_t numSlices, double Xdropcutoff, const std::vector<std::tuple<size_t, int, int>>& forbiddenNodes, AlignerGraphsizedState& reusableState) const
	{
		assert(initialSlice.j == (size_t)-WordConfiguration<Word>::WordSize);
		assert(initialSlice.j + numSlices * WordConfiguration<Word>::WordSize <= sequence.size() + WordConfiguration<Word>::WordSize);
		DPTable result;
		result.slices.reserve(numSlices + 1);
		size_t cellsProcessed = 0;
		std::vector<size_t> partOfComponent;
		{
			for (auto node : initialSlice.scores)
			{
				reusableState.previousBand[node.first] = true;
			}
		}
		DPSlice lastSlice = initialSlice;
		result.slices.push_back(initialSlice);
		ScoreType bestXScore = initialSlice.maxExactEndposScore;
		assert(bestXScore != std::numeric_limits<ScoreType>::min());
		std::vector<ProcessedSeedHit> fakeSeeds;
		WordSlice fakeSlice { WordConfiguration<Word>::AllZeros, WordConfiguration<Word>::AllZeros, std::numeric_limits<ScoreType>::max() };
		auto nodeAllowanceVectors = getNodeAllowanceVectors(forbiddenNodes, numSlices);
		for (size_t slice = 0; slice < numSlices; slice++)
		{
			int bandwidth = params.alignmentBandwidth;
#ifdef SLICEVERBOSE
			auto timeStart = std::chrono::system_clock::now();
#endif
			DPSlice newSlice;
			fixAllowedNodes(reusableState.allowedBigraphNodesThisSlice, nodeAllowanceVectors, slice);
			if (reusableState.componentQueue.valid())
			{
				newSlice = pickMethodAndExtendFill(sequence, lastSlice, reusableState.previousBand, reusableState.currentBand, reusableState.componentQueue, bandwidth, fakeSeeds, std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), fakeSlice, reusableState.hasSeedStart, false, reusableState.allowedBigraphNodesThisSlice);
			}
			else
			{
				newSlice = pickMethodAndExtendFill(sequence, lastSlice, reusableState.previousBand, reusableState.currentBand, reusableState.calculableQueue, bandwidth, fakeSeeds, std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), fakeSlice, reusableState.hasSeedStart, false, reusableState.allowedBigraphNodesThisSlice);
			}
#ifdef SLICEVERBOSE
			auto timeEnd = std::chrono::system_clock::now();
			auto time = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart).count();
			std::cerr << "slice " << slice << " bandwidth " << bandwidth << " minscore " << newSlice.minScore << " diff " << (newSlice.minScore - lastSlice.minScore) << " time " << time << " nodes " << newSlice.scores.size() << " slices " << newSlice.numCells << " nodesprocessed " << newSlice.nodesProcessed << " cellsprocessed " << newSlice.cellsProcessed << " overhead " << (100 * (int)(newSlice.cellsProcessed - newSlice.numCells) / (int)(newSlice.numCells)) << "%";
#endif
			assert(newSlice.minScore != std::numeric_limits<ScoreType>::max());
			assert(newSlice.minScoreNode != std::numeric_limits<LengthType>::max());
			assert(newSlice.minScoreNodeOffset != std::numeric_limits<LengthType>::max());
			assert(newSlice.scores.hasNode(newSlice.minScoreNode));
			assert(newSlice.minScoreNodeOffset < params.graph.NodeLength(newSlice.minScoreNode));
			assert(newSlice.maxExactEndposScore != std::numeric_limits<ScoreType>::min());

			if (newSlice.maxExactEndposScore > bestXScore) bestXScore = newSlice.maxExactEndposScore;

			assert(newSlice.j == lastSlice.j + WordConfiguration<Word>::WordSize);

			cellsProcessed += newSlice.cellsProcessed;

			if (newSlice.cellsProcessed >= params.maxCellsPerSlice)
			{
				newSlice.scoresNotValid = true;
			}

			if (newSlice.maxExactEndposScore < bestXScore - Xdropcutoff)
			{
				for (auto node : lastSlice.scores)
				{
					assert(reusableState.previousBand[node.first]);
					reusableState.previousBand[node.first] = false;
				}
				for (auto node : newSlice.scores)
				{
					assert(reusableState.currentBand[node.first]);
					reusableState.currentBand[node.first] = false;
				}
				lastSlice.scoresVectorMap.removeVectorArray();
				newSlice.scoresVectorMap.removeVectorArray();
				break;
			}

#ifdef SLICEVERBOSE
			std::cerr << std::endl;
#endif

			result.slices.push_back(newSlice.getMapSlice());
			for (auto node : lastSlice.scores)
			{
				assert(reusableState.previousBand[node.first]);
				reusableState.previousBand[node.first] = false;
			}
			assert(newSlice.minScore != std::numeric_limits<ScoreType>::max());
			assert(newSlice.minScore >= lastSlice.minScore);
			if (slice == numSlices - 1)
			{
				for (auto node : newSlice.scores)
				{
					assert(reusableState.currentBand[node.first]);
					reusableState.currentBand[node.first] = false;
				}
			}
			else
			{
				std::swap(reusableState.previousBand, reusableState.currentBand);
			}
			lastSlice.scoresVectorMap.removeVectorArray();
			lastSlice = std::move(newSlice);
		}
		lastSlice.scoresVectorMap.removeVectorArray();

		assert(result.slices.size() <= numSlices + 1);

#ifdef EXTRACORRECTNESSASSERTIONS
		assert(reusableState.calculableQueue.size() == 0);
		for (size_t i = 0; i < reusableState.currentBand.size(); i++)
		{
			assert(!reusableState.currentBand[i]);
			assert(!reusableState.previousBand[i]);
		}
#endif

#ifndef NDEBUG
		if (result.slices.size() > 0)
		{
			for (size_t i = 1; i < result.slices.size(); i++)
			{
				assert(result.slices[i].j == result.slices[i-1].j + WordConfiguration<Word>::WordSize);
			}
			for (size_t i = 1; i < result.slices.size(); i++)
			{
				assert(result.slices[i].minScore >= result.slices[i-1].minScore);
			}
		}
#endif
		return result;
	}

	std::pair<std::vector<std::vector<size_t>>, std::vector<std::vector<size_t>>> getNodeAllowanceVectors(const std::vector<std::tuple<size_t, int, int>>& forbiddenSpans, size_t numSlices) const
	{
		std::vector<std::vector<size_t>> newlyAllowed;
		std::vector<std::vector<size_t>> newlyForbidden;
		newlyForbidden.resize(numSlices);
		newlyAllowed.resize(numSlices);
		for (auto t : forbiddenSpans)
		{
			size_t forbidslice = 0;
			if (std::get<1>(t) >= 1)
			{
				forbidslice = (std::get<1>(t) + 63) / 64;
				if (forbidslice >= numSlices) forbidslice = numSlices-1;
			}
			size_t allowslice = 0;
			if (std::get<2>(t) >= 1)
			{
				allowslice = std::get<2>(t) / 64;
				if (allowslice >= numSlices) allowslice = numSlices-1;
			}
			if (forbidslice == allowslice) continue;
			assert(allowslice > forbidslice);
			assert(forbidslice < newlyForbidden.size());
			assert(allowslice < newlyAllowed.size());
			newlyForbidden[forbidslice].emplace_back(std::get<0>(t));
			newlyAllowed[allowslice].emplace_back(std::get<0>(t));
		}
		return std::make_pair(std::move(newlyAllowed), std::move(newlyForbidden));
	}

	void fixAllowedNodes(std::vector<bool>& allowedBigraphNodesThisSlice, const std::pair<std::vector<std::vector<size_t>>, std::vector<std::vector<size_t>>>& nodeAllowanceVectors, const size_t slice) const
	{
		assert(slice < nodeAllowanceVectors.first.size());
		assert(slice < nodeAllowanceVectors.second.size());
		for (size_t node : nodeAllowanceVectors.first[slice])
		{
			assert(node < allowedBigraphNodesThisSlice.size());
			allowedBigraphNodesThisSlice[node] = true;
		}
		for (size_t node : nodeAllowanceVectors.second[slice])
		{
			assert(node < allowedBigraphNodesThisSlice.size());
			allowedBigraphNodesThisSlice[node] = false;
		}
	}

  /*
   * zkn
   * This function aligns an entire read. It takes as an input the entire
   * cluster, and extends each seed in the cluster individually, starting with
   * the best candidate seeds. Seeds which are overlapped by previously computed
   * alignments are dropped to avoid redundant work
   */
	DPTable getMultiseedSlices(
      const std::string_view& sequence,
      const DPSlice& initialSlice,
      size_t numSlices,
      AlignerGraphsizedState& reusableState,
      const std::vector<ProcessedSeedHit>& seedHits,
      const std::vector<ScoreType>& sliceMaxScores,
      const std::vector<std::tuple<size_t, int, int>>& forbiddenNodes) const
	{

		assert(reusableState.componentQueue.valid());
		assert(initialSlice.j == (size_t)-WordConfiguration<Word>::WordSize);
		assert(initialSlice.j + numSlices * WordConfiguration<Word>::WordSize <= sequence.size() + WordConfiguration<Word>::WordSize);
		assert(numSlices >= 1);
		assert(numSlices == 1 || initialSlice.j + (numSlices-1) * WordConfiguration<Word>::WordSize <= sequence.size());
		DPTable result;
		result.slices.reserve(numSlices + 1);
		size_t cellsProcessed = 0;
		DPSlice lastSlice = initialSlice;
		result.slices.push_back(initialSlice);
		size_t lastSeedHit = 0;
		ScoreType XDropCurrentBest = 0;
		assert(params.Xdropcutoff > 0);
		auto nodeAllowanceVectors = getNodeAllowanceVectors(forbiddenNodes, numSlices);
		for (size_t slice = 0; slice < numSlices; slice++)
		{
			int bandwidth = params.alignmentBandwidth;
#ifdef SLICEVERBOSE
			auto timeStart = std::chrono::system_clock::now();
#endif
			fixAllowedNodes(reusableState.allowedBigraphNodesThisSlice, nodeAllowanceVectors, slice);
			size_t nextSeedHit = lastSeedHit;
			assert(nextSeedHit == seedHits.size() || seedHits[nextSeedHit].seqPos / WordConfiguration<Word>::WordSize >= (lastSlice.j + WordConfiguration<Word>::WordSize) / WordConfiguration<Word>::WordSize);
			while (nextSeedHit < seedHits.size() && seedHits[nextSeedHit].seqPos / WordConfiguration<Word>::WordSize == (lastSlice.j + WordConfiguration<Word>::WordSize) / WordConfiguration<Word>::WordSize)
			{
				nextSeedHit += 1;
			}
			assert(nextSeedHit == seedHits.size() || seedHits[nextSeedHit].seqPos / WordConfiguration<Word>::WordSize > (lastSlice.j + WordConfiguration<Word>::WordSize) / WordConfiguration<Word>::WordSize);
			size_t seqOffset = lastSlice.j + WordConfiguration<Word>::WordSize;
			WordSlice seedSlice = BV::getSeedSlice(seqOffset, sequence.size(), params);
			assert(seedSlice.maxXScore(seqOffset, params.XscoreErrorCost) >= -(ScoreType)WordConfiguration<Word>::WordSize*100);
			assert(seedSlice.maxXScore(seqOffset, params.XscoreErrorCost) <= (ScoreType)WordConfiguration<Word>::WordSize*100);
			ScoreType possibleScoreRemaining = sequence.size() - seqOffset;
			assert(possibleScoreRemaining > 0);
			// // can't get an alignment which would be backtraced. fake set the seed set to be empty
			// if (possibleScoreRemaining < sliceMaxScores[slice]) lastSeedHit = nextSeedHit;
			DPSlice newSlice = pickMethodAndExtendFill(sequence, lastSlice, reusableState.previousBand, reusableState.currentBand, reusableState.componentQueue, bandwidth, seedHits, lastSeedHit, nextSeedHit, seedSlice, reusableState.hasSeedStart, true, reusableState.allowedBigraphNodesThisSlice);
			lastSeedHit = nextSeedHit;
			XDropCurrentBest = std::max(XDropCurrentBest, newSlice.maxExactEndposScore);
			if (newSlice.maxExactEndposScore < XDropCurrentBest - params.Xdropcutoff)
			{
				for (auto node : newSlice.scores)
				{
					reusableState.currentBand[node.first] = false;
				}
				newSlice.minScore = seedSlice.getScoreBeforeStart();
				newSlice.minScoreNode = std::numeric_limits<LengthType>::max();
				newSlice.maxExactEndposScore = std::numeric_limits<ScoreType>::min();
				newSlice.maxExactEndposNode = std::numeric_limits<LengthType>::max();
				newSlice.nodeMaxExactEndposScore.clear();
				newSlice.seedstartNodes.clear();
				newSlice.scores.clear();
				XDropCurrentBest = 0;
			}
#ifdef SLICEVERBOSE
			auto timeEnd = std::chrono::system_clock::now();
			auto time = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart).count();
			std::cerr << "slice " << slice << " bandwidth " << bandwidth << " minscore " << newSlice.minScore << " diff " << (newSlice.minScore - lastSlice.minScore) << " maxX " << newSlice.maxExactEndposScore << " time " << time << " nodes " << newSlice.scores.size() << " slices " << newSlice.numCells << " nodesprocessed " << newSlice.nodesProcessed << " cellsprocessed " << newSlice.cellsProcessed << " overhead " << (100 * (int)(newSlice.cellsProcessed - newSlice.numCells) / (std::max((int)newSlice.numCells, (int)1))) << "%";
#endif

			assert(newSlice.j == lastSlice.j + WordConfiguration<Word>::WordSize);

			cellsProcessed += newSlice.cellsProcessed;

			if (newSlice.cellsProcessed >= params.maxCellsPerSlice)
			{
				newSlice.scoresNotValid = true;
			}

#ifdef SLICEVERBOSE
			std::cerr << std::endl;
#endif

			result.slices.push_back(newSlice.getMapSlice());
			for (auto node : lastSlice.scores)
			{
				assert(reusableState.previousBand[node.first]);
				reusableState.previousBand[node.first] = false;
			}
			if (slice == numSlices - 1)
			{
				for (auto node : newSlice.scores)
				{
					assert(reusableState.currentBand[node.first]);
					reusableState.currentBand[node.first] = false;
				}
			}
			else
			{
				std::swap(reusableState.previousBand, reusableState.currentBand);
			}
			lastSlice.scoresVectorMap.removeVectorArray();
			lastSlice = std::move(newSlice);
		}
		assert(lastSeedHit == seedHits.size());
		lastSlice.scoresVectorMap.removeVectorArray();

		assert(result.slices.size() == numSlices + 1);

#ifdef EXTRACORRECTNESSASSERTIONS
		assert(reusableState.calculableQueue.size() == 0);
		for (size_t i = 0; i < reusableState.currentBand.size(); i++)
		{
			assert(!reusableState.currentBand[i]);
			assert(!reusableState.previousBand[i]);
		}
#endif

#ifndef NDEBUG
		if (result.slices.size() > 0)
		{
			for (size_t i = 1; i < result.slices.size(); i++)
			{
				assert(result.slices[i].j == result.slices[i-1].j + WordConfiguration<Word>::WordSize);
			}
		}
#endif
		return result;
	}

};

#endif
