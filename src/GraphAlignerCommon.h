#ifndef GraphAlignerCommon_h
#define GraphAlignerCommon_h

#include <vector>
#include "AlignmentGraph.h"
#include "ArrayPriorityQueue.h"
#include "ComponentPriorityQueue.h"
#include "NodeSlice.h"
#include "WordSlice.h"
#include "vg.pb.h"

#include "nlohmann/json.hpp"

using nJson = nlohmann::json;

template <typename LengthType, typename ScoreType, typename Word>
class GraphAlignerCommon
{
public:
	class NodeWithPriority
	{
	public:
		NodeWithPriority(LengthType node, size_t offset, size_t endOffset, ScoreType priority) : node(node), offset(offset), endOffset(endOffset), priority(priority) {}
		bool operator>(const NodeWithPriority& other) const
		{
		return priority > other.priority;
		}
		bool operator<(const NodeWithPriority& other) const
		{
		return priority < other.priority;
		}
		LengthType node;
		size_t offset;
		size_t endOffset;
		ScoreType priority;
	};
	class EdgeWithPriority
	{
	public:
		EdgeWithPriority(LengthType target, ScoreType priority, WordSlice<LengthType, ScoreType, Word> incoming, bool skipFirst) : target(target), priority(priority), incoming(incoming), skipFirst(skipFirst), slice(0), forceCalculation(false) {}
		EdgeWithPriority(LengthType target, ScoreType priority, WordSlice<LengthType, ScoreType, Word> incoming, bool skipFirst, size_t slice) : target(target), priority(priority), incoming(incoming), skipFirst(skipFirst), slice(slice), forceCalculation(false) {}
		bool operator>(const EdgeWithPriority& other) const
		{
			return priority > other.priority;
		}
		bool operator<(const EdgeWithPriority& other) const
		{
			return priority < other.priority;
		}
		LengthType target;
		ScoreType priority;
		WordSlice<LengthType, ScoreType, Word> incoming;
		bool skipFirst;
		size_t slice;
		bool forceCalculation;
	};
	class AlignerGraphsizedState
	{
	public:
		AlignerGraphsizedState(const AlignmentGraph& graph, size_t maxBandwidth) :
		componentQueue(),
		calculableQueue(),
		currentBand(),
		previousBand(),
		hasSeedStart(),
		allowedBigraphNodesThisSlice(),
		bigraphNodeForbiddenSpans()
		{
			componentQueue.initialize(graph.ComponentSize());
			calculableQueue.initialize(WordConfiguration<Word>::WordSize * (WordConfiguration<Word>::WordSize + maxBandwidth + 1) + maxBandwidth + 1, graph.NodeSize());
			currentBand.resize(graph.NodeSize(), false);
			previousBand.resize(graph.NodeSize(), false);
			hasSeedStart.resize(graph.NodeSize(), false);
			allowedBigraphNodesThisSlice.resize(graph.BigraphNodeCount(), true);
		}
    //zkn It appears this pretty much sets it back to the way it was when it was
    //first constructed, which means it's safe to use in the kernel for loop
    //before each read. I believe this is usually done outside the kernel, which
    //is why it must be called manually before
		void clear()
		{
			componentQueue.clear();
			calculableQueue.clear();
			currentBand.assign(currentBand.size(), false);
			previousBand.assign(previousBand.size(), false);
			hasSeedStart.assign(hasSeedStart.size(), false);
			allowedBigraphNodesThisSlice.assign(allowedBigraphNodesThisSlice.size(), true);
			bigraphNodeForbiddenSpans.clear();
		}
		ComponentPriorityQueue<EdgeWithPriority, true> componentQueue;
		ArrayPriorityQueue<EdgeWithPriority, true> calculableQueue;
		std::vector<bool> currentBand;
		std::vector<bool> previousBand;
		std::vector<bool> hasSeedStart;
		std::vector<bool> allowedBigraphNodesThisSlice;
		std::vector<std::tuple<size_t, int, int>> bigraphNodeForbiddenSpans;
	};
	using MatrixPosition = AlignmentGraph::MatrixPosition;

  //fwd declaration
  class SerializableParams;
	class Params
	{
	public:
		Params(LengthType alignmentBandwidth, const AlignmentGraph& graph, size_t maxCellsPerSlice, bool quietMode, double preciseClippingIdentityCutoff, ScoreType Xdropcutoff, double multimapScoreFraction, ScoreType clipAmbiguousEnds, size_t maxTraceCount) :
		alignmentBandwidth(alignmentBandwidth),
		graph(graph),
		maxCellsPerSlice(maxCellsPerSlice),
		quietMode(quietMode),
		XscoreErrorCost(100 * (preciseClippingIdentityCutoff / (1.0 - preciseClippingIdentityCutoff) + 1.0)),
		Xdropcutoff(Xdropcutoff),
		multimapScoreFraction(multimapScoreFraction),
		clipAmbiguousEnds(clipAmbiguousEnds),
		discardCigar(false),
		maxTraceCount(maxTraceCount)
		{
		}
    // Constructor that takes a SerializableParams object
    Params(const SerializableParams& params) :
        alignmentBandwidth(params.alignmentBandwidth),
        graph(params.graph),
        maxCellsPerSlice(params.maxCellsPerSlice),
        quietMode(params.quietMode),
        XscoreErrorCost(params.XscoreErrorCost),
        Xdropcutoff(params.Xdropcutoff),
        multimapScoreFraction(params.multimapScoreFraction),
        clipAmbiguousEnds(params.clipAmbiguousEnds),
        maxTraceCount(params.maxTraceCount),
        discardCigar(params.discardCigar)
    {
    }
		const LengthType alignmentBandwidth;
    //zkn This is the only attribute that requires more complicated pickling
		const AlignmentGraph& graph;
		const size_t maxCellsPerSlice;
		const bool quietMode;
		const ScoreType XscoreErrorCost;
		const ScoreType Xdropcutoff;
		const double multimapScoreFraction;
		const ScoreType clipAmbiguousEnds;
		const size_t maxTraceCount;
		bool discardCigar;
	};

  class SerializableParams {
  private:
    /* zkn boost serialization stuff */
    friend class boost::serialization::access;
    /*
     * Not actually serializing any members, just want the template params
     */
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version){
		  ar & alignmentBandwidth;
      ar & graph;
      ar & maxCellsPerSlice;
      ar & quietMode;
      ar & XscoreErrorCost;
      ar & Xdropcutoff;
      ar & multimapScoreFraction;
      ar & clipAmbiguousEnds;
      ar & maxTraceCount;
      ar & discardCigar;
    }
  public:
      // Constructor that takes a Params object
      SerializableParams(const Params& params) :
          alignmentBandwidth(params.alignmentBandwidth),
          graph(params.graph),
          maxCellsPerSlice(params.maxCellsPerSlice),
          quietMode(params.quietMode),
          XscoreErrorCost(params.XscoreErrorCost),
          Xdropcutoff(params.Xdropcutoff),
          multimapScoreFraction(params.multimapScoreFraction),
          clipAmbiguousEnds(params.clipAmbiguousEnds),
          discardCigar(params.discardCigar),
          maxTraceCount(params.maxTraceCount)
      {
      }
      //default
      SerializableParams() :
          alignmentBandwidth(0),      
          graph(AlignmentGraph()),    
          maxCellsPerSlice(0),        
          quietMode(false),           
          XscoreErrorCost(0),         
          Xdropcutoff(0),             
          multimapScoreFraction(0.0), 
          clipAmbiguousEnds(0),       
          maxTraceCount(0),           
          discardCigar(false)         
      {
      }

      // Members without const
      LengthType alignmentBandwidth;
      AlignmentGraph graph;
      size_t maxCellsPerSlice;
      bool quietMode;
      ScoreType XscoreErrorCost;
      ScoreType Xdropcutoff;
      double multimapScoreFraction;
      ScoreType clipAmbiguousEnds;
      size_t maxTraceCount;
      bool discardCigar;
  };


	struct TraceItem
	{
		TraceItem() :
		DPposition(),
		nodeSwitch(false),
		sequenceCharacter('-'),
		graphCharacter('-')
		{}
		TraceItem(MatrixPosition DPposition, bool nodeSwitch, char sequenceCharacter, char graphCharacter) :
		DPposition(DPposition),
		nodeSwitch(nodeSwitch),
		sequenceCharacter(sequenceCharacter),
		graphCharacter(graphCharacter)
		{}
		TraceItem(MatrixPosition DPposition, bool nodeSwitch, const std::string& seq, const AlignmentGraph& graph) :
		DPposition(DPposition),
		nodeSwitch(nodeSwitch),
		sequenceCharacter(DPposition.seqPos < seq.size() ? seq[DPposition.seqPos] : '-'),
		graphCharacter(graph.NodeSequences(DPposition.node, DPposition.nodeOffset))
		{}
		TraceItem(MatrixPosition DPposition, bool nodeSwitch, const std::string_view& seq, const AlignmentGraph& graph) :
		DPposition(DPposition),
		nodeSwitch(nodeSwitch),
		sequenceCharacter(DPposition.seqPos < seq.size() ? seq[DPposition.seqPos] : '-'),
		graphCharacter(graph.NodeSequences(DPposition.node, DPposition.nodeOffset))
		{}
		bool operator==(const TraceItem& other) const
		{
			return DPposition == other.DPposition && nodeSwitch == other.nodeSwitch && sequenceCharacter == other.sequenceCharacter && graphCharacter == other.graphCharacter;
		}
		MatrixPosition DPposition;
		bool nodeSwitch;
		char sequenceCharacter;
		char graphCharacter;
	};
	class OnewayTrace
	{
	public:
		OnewayTrace() :
		trace(),
		score(0)
		{
		}
		// force move semantics because copying is very slow and unnecessary
		OnewayTrace(const OnewayTrace& other) = delete;
		OnewayTrace(OnewayTrace&& other) = default;
		OnewayTrace& operator=(const OnewayTrace& other) = delete;
		OnewayTrace& operator=(OnewayTrace&& other) = default;
		static OnewayTrace TraceFailed()
		{
			OnewayTrace result;
			result.score = std::numeric_limits<ScoreType>::max();
			return result;
		}
    /*
     * zkn
     * to string method for oneway trace. This is used withing GraphAligner.h in
     * order to dump the output to a file. We use json to represent the string,
     * but it is not a complete serialization of the file. Niether does it
     * conform to any standardized genomics file formats (sam, gaf, ...).
     * We include the following data structures.
     * {
     *   "score"      : alignment score,
     *   "seqChars"   : sequence chars,
     *   "nodeId"     : node that the char matched to,
     *   "nodeOffset" : offset of the node it matched to,
     *   "graphChars" : corresponding characters in the graph
     * }
     */
    nJson toJson(){
      nJson o;
      o["score"]      = score;
      o["seqChars"]   = nJson::array();
      o["nodeId"]     = nJson::array();
      o["nodeOffset"]     = nJson::array();
      o["graphChars"] = nJson::array();
      for(TraceItem& t : trace){
        o["seqChars"].push_back(t.sequenceCharacter);
        o["nodeId"].push_back(t.DPposition.node);
        o["nodeOffset"].push_back(t.DPposition.nodeOffset);
        o["graphChars"].push_back(t.graphCharacter);
      }
      return o;
    }
		bool failed() const
		{
			return score == std::numeric_limits<ScoreType>::max();
		}
		ScoreType alignmentXScore(ScoreType XscoreErrorCost) const
		{
			if (failed()) return 0;
			assert(trace.size() > 0);
			assert(trace[0].DPposition.seqPos >= trace.back().DPposition.seqPos);
			return (ScoreType)(trace[0].DPposition.seqPos - trace.back().DPposition.seqPos + 1)*100 - XscoreErrorCost * (ScoreType)score;
		}
		std::vector<TraceItem> trace;
		ScoreType score;
	};
	class Trace
	{
	public:
		OnewayTrace forward;
		OnewayTrace backward;
	};
#ifdef NDEBUG
	__attribute__((always_inline))
#endif
	static bool characterMatch(char sequenceCharacter, char graphCharacter)
	{
		if (sequenceCharacter == graphCharacter) return true;
		switch(sequenceCharacter)
		{
			case 'a':
			case 'A':
				return ambiguousMatch(graphCharacter, 'A');
			case 'c':
			case 'C':
				return ambiguousMatch(graphCharacter, 'C');
			case 'g':
			case 'G':
				return ambiguousMatch(graphCharacter, 'G');
			case 't':
			case 'T':
				return ambiguousMatch(graphCharacter, 'T');
			case '-':
				return false;
		}
		return (ambiguousMatch(sequenceCharacter, 'A') && ambiguousMatch(graphCharacter, 'A'))
		|| (ambiguousMatch(sequenceCharacter, 'C') && ambiguousMatch(graphCharacter, 'C'))
		|| (ambiguousMatch(sequenceCharacter, 'G') && ambiguousMatch(graphCharacter, 'G'))
		|| (ambiguousMatch(sequenceCharacter, 'T') && ambiguousMatch(graphCharacter, 'T'));
	}
#ifdef NDEBUG
	__attribute__((always_inline))
#endif
	static bool ambiguousMatch(char ambiguousChar, char exactChar)
	{
		assert(exactChar == 'A' || exactChar == 'T' || exactChar == 'C' || exactChar == 'G');
		switch(ambiguousChar)
		{
			case '-':
				return false;
			case 'A':
			case 'a':
				return exactChar == 'A';
			break;
			case 'u':
			case 'U':
			case 'T':
			case 't':
				return exactChar == 'T';
			break;
			case 'C':
			case 'c':
				return exactChar == 'C';
			break;
			case 'G':
			case 'g':
				return exactChar == 'G';
			break;
			case 'N':
			case 'n':
				return true;
			break;
			case 'R':
			case 'r':
				return exactChar == 'A' || exactChar == 'G';
			break;
			case 'Y':
			case 'y':
				return exactChar == 'C' || exactChar == 'T';
			break;
			case 'K':
			case 'k':
				return exactChar == 'G' || exactChar == 'T';
			break;
			case 'M':
			case 'm':
				return exactChar == 'C' || exactChar == 'A';
			break;
			case 'S':
			case 's':
				return exactChar == 'C' || exactChar == 'G';
			break;
			case 'W':
			case 'w':
				return exactChar == 'A' || exactChar == 'T';
			break;
			case 'B':
			case 'b':
				return exactChar == 'C' || exactChar == 'G' || exactChar == 'T';
			break;
			case 'D':
			case 'd':
				return exactChar == 'A' || exactChar == 'G' || exactChar == 'T';
			break;
			case 'H':
			case 'h':
				return exactChar == 'A' || exactChar == 'C' || exactChar == 'T';
			break;
			case 'V':
			case 'v':
				return exactChar == 'A' || exactChar == 'C' || exactChar == 'G';
			break;
			default:
				assert(false);
				std::abort();
				return false;
			break;
		}
	}
};

class AlignmentResult
{
public:
	AlignmentResult() :
		alignments(),
		seedsExtended(0)
	{}
	class AlignmentItem
	{
	public:
		AlignmentItem() :
		corrected(),
		alignment(),
		trace(),
		seedGoodness(0),
		cellsProcessed(0),
		elapsedMilliseconds(0),
		alignmentStart(0),
		alignmentEnd(0),
		alignmentScore(std::numeric_limits<size_t>::max()),
		alignmentXScore(-1),
		mappingQuality(255)
		{}
		AlignmentItem(GraphAlignerCommon<size_t, int64_t, uint64_t>::OnewayTrace&& trace, size_t cellsProcessed, size_t ms) :
		corrected(),
		alignment(),
		trace(),
		cellsProcessed(cellsProcessed),
		elapsedMilliseconds(ms),
		alignmentStart(0),
		alignmentEnd(0),
		alignmentScore(std::numeric_limits<size_t>::max()),
		alignmentXScore(-1),
		mappingQuality(255)
		{
			this->trace = std::make_shared<GraphAlignerCommon<size_t, int64_t, uint64_t>::OnewayTrace>();
			*this->trace = std::move(trace);
		}
		bool alignmentFailed() const
		{
			return alignmentEnd == alignmentStart;
		}
		size_t alignmentLength() const
		{
			return alignmentEnd - alignmentStart;
		}
		std::string corrected;
		std::string GAFline;
		std::shared_ptr<vg::Alignment> alignment;
		std::shared_ptr<GraphAlignerCommon<size_t, int64_t, uint64_t>::OnewayTrace> trace;
		size_t seedGoodness;
		size_t cellsProcessed;
		size_t elapsedMilliseconds;
		size_t alignmentStart;
		size_t alignmentEnd;
		size_t alignmentScore;
		double alignmentXScore;
		int mappingQuality;
	};
	std::vector<AlignmentItem> alignments;
	size_t seedsExtended;
	std::string readName;
};

#endif
