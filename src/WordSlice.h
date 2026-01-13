#ifndef WordSlice_h
#define WordSlice_h

#include "DebugFlags.h"

template <typename LengthType, typename ScoreType, typename Word>
class GraphAlignerBitvectorCommon;

template <typename Word>
class WordConfiguration
{
};

template <>
class WordConfiguration<uint64_t>
{
public:
	static constexpr int WordSize = 64;
	//number of bits per chunk
	//prefix sum differences are calculated in chunks of log w bits
	static constexpr int ChunkBits = 8;
	static constexpr uint64_t AllZeros = 0x0000000000000000;
	static constexpr uint64_t AllOnes = 0xFFFFFFFFFFFFFFFF;
	static constexpr uint64_t LastBit = 0x8000000000000000;
	//positions of the sign bits for each chunk
	static constexpr uint64_t SignMask = 0x8080808080808080;
	//constant for multiplying the chunk popcounts into prefix sums
	//this should be 1 at the start of each chunk
	static constexpr uint64_t PrefixSumMultiplierConstant = 0x0101010101010101;
	//positions of the least significant bits for each chunk
	static constexpr uint64_t LSBMask = 0x0101010101010101;

#ifdef NOBUILTINPOPCOUNT
	static int popcount(uint64_t x)
	{
		//https://en.wikipedia.org/wiki/Hamming_weight
		x -= (x >> 1) & 0x5555555555555555;
		x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333);
		x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0f;
		return (x * 0x0101010101010101) >> 56;
	}
#else
	static int popcount(uint64_t x)
	{
		//https://gcc.gnu.org/onlinedocs/gcc-4.8.4/gcc/X86-Built-in-Functions.html
		// return __builtin_popcountll(x);
		//for some reason __builtin_popcount takes 21 instructions so call assembly directly
		__asm__("popcnt %0, %0" : "+r" (x));
		return x;
	}
#endif


	static uint64_t ChunkPopcounts(uint64_t value)
	{
		uint64_t x = value;
		x -= (x >> 1) & 0x5555555555555555;
		x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333);
		x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0f;
		return x;
	}

	static int BitPosition(uint64_t low, uint64_t high, int rank)
	{
		assert(rank >= 0);
		if (popcount(low) <= rank)
		{
			rank -= popcount(low);
			return 64 + BitPosition(high, rank);
		}
		return BitPosition(low, rank);
	}

	static int BitPosition(uint64_t number, int rank)
	{
		// // https://stackoverflow.com/questions/45482787/how-to-efficiently-find-the-n-th-set-bit
		// uint64_t j = _pdep_u64(1 << rank, number);
		// return (__builtin_ctz(j));
		uint64_t bytes = ChunkPopcounts(number);
		//cumulative popcount of each byte
		uint64_t cumulative = bytes * PrefixSumMultiplierConstant;
		//spread the rank into each byte
		uint64_t rankFinder = ((rank + 1) & 0xFF) * PrefixSumMultiplierConstant;
		//rankMask's msb will be 0 if the c. popcount at that byte is < rank, or 1 if >= rank
		uint64_t rankMask = (cumulative | SignMask) - rankFinder;
		//the total number of ones in rankMask is the number of bytes whose c. popcount is >= rank
		//8 - that is the number of bytes whose c. popcount is < rank
		int smallerBytes = 8 - ((((rankMask & SignMask) >> 7) * PrefixSumMultiplierConstant) >> 56);
		assert(smallerBytes < 8);
		//the bit position will be inside this byte
		uint64_t interestingByte = (number >> (smallerBytes * 8)) & 0xFF;
		if (smallerBytes > 0) rank -= (cumulative >> ((smallerBytes - 1) * 8)) & 0xFF;
		assert(rank >= 0 && rank < 8);
		//spread the 1's from interesting byte to each byte
		//first put every pair of bits into each 2-byte boundary
		//then select only those pairs
		//then spread the pairs into each byte boundary
		//and select the ones
		uint64_t spreadBits = (((interestingByte * 0x0000040010004001) & 0x0003000300030003) * 0x0000000000000081) & 0x0101010101010101;
/*
0000 0000  0000 0000  0000 0000  0000 0000  0000 0000  0000 0000  0000 0000  abcd efgh
0000 0000  0000 00ab  cdef gh00  0000 abcd  efgh 0000  00ab cdef  gh00 0000  abcd efgh  * 0x0000040010004001
0000 0000  0000 00ab  0000 0000  0000 00cd  0000 0000  0000 00ef  0000 0000  0000 00gh  & 0x0003000300030003
0000 000a  b000 00ab  0000 000c  d000 00cd  0000 000e  f000 00ef  0000 000g  h000 00gh  * 0x0000000000000081
0000 000a  0000 000b  0000 000c  0000 000d  0000 000e  0000 000f  0000 000g  0000 000h  & 0x0101010101010101
*/
		//find the position from the bits the same way as from the bytes
		uint64_t cumulativeBits = spreadBits * PrefixSumMultiplierConstant;
		uint64_t bitRankFinder = ((rank + 1) & 0xFF) * PrefixSumMultiplierConstant;
		uint64_t bitRankMask = (cumulativeBits | SignMask) - bitRankFinder;
		int smallerBits = 8 - ((((bitRankMask & SignMask) >> 7) * PrefixSumMultiplierConstant) >> 56);
		assert(smallerBits >= 0);
		assert(smallerBits < 8);
		return smallerBytes * 8 + smallerBits;
	}

	static uint64_t MortonHigh(uint64_t left, uint64_t right)
	{
		return Interleave(left >> 32, right >> 32);
	}

	static uint64_t MortonLow(uint64_t left, uint64_t right)
	{
		return Interleave(left & 0xFFFFFFFF, right & 0xFFFFFFFF);
	}

	//http://graphics.stanford.edu/~seander/bithacks.html#InterleaveBMN
	static uint64_t Interleave(uint64_t x, uint64_t y)
	{
		assert(x == (x & 0xFFFFFFFF));
		assert(y == (y & 0xFFFFFFFF));
		static const uint64_t B[] = {0x5555555555555555, 0x3333333333333333, 0x0F0F0F0F0F0F0F0F, 0x00FF00FF00FF00FF, 0x0000FFFF0000FFFF};
		static const uint64_t S[] = {1, 2, 4, 8, 16};

		x = (x | (x << S[4])) & B[4];
		x = (x | (x << S[3])) & B[3];
		x = (x | (x << S[2])) & B[2];
		x = (x | (x << S[1])) & B[1];
		x = (x | (x << S[0])) & B[0];

		y = (y | (y << S[4])) & B[4];
		y = (y | (y << S[3])) & B[3];
		y = (y | (y << S[2])) & B[2];
		y = (y | (y << S[1])) & B[1];
		y = (y | (y << S[0])) & B[0];

		return x | (y << 1);
	}
};

//uncomment if there's an undefined reference with -O0. why?
// constexpr uint64_t WordConfiguration<uint64_t>::AllZeros;
// constexpr uint64_t WordConfiguration<uint64_t>::AllOnes;

template <typename LengthType, typename ScoreType, typename Word>
class WordSlice
{
public:
	WordSlice() :
	VP(0),
	VN(0),
	scoreEnd(0)
	{}
	WordSlice(Word VP, Word VN, ScoreType scoreEnd) :
	VP(VP),
	VN(VN),
	scoreEnd(scoreEnd)
	{}
	Word VP;
	Word VN;
	ScoreType scoreEnd;

#ifdef NDEBUG
	__attribute__((always_inline))
#endif
	WordSlice mergeWith(const WordSlice& other) const
	{
		auto result = mergeTwoSlices(*this, other);
		return result;
	}

	ScoreType getValue(int row) const
	{
		Word mask = WordConfiguration<Word>::AllZeros;
		if (row < WordConfiguration<Word>::WordSize-1) mask = (WordConfiguration<Word>::AllOnes << (row + 1));
		auto value = scoreEnd + WordConfiguration<Word>::popcount(VN & mask) - WordConfiguration<Word>::popcount(VP & mask);
#ifdef EXTRACORRECTNESSASSERTIONS
		assert(value == getValueStartBased(row));
#endif
		return value;
	}

	ScoreType getValueStartBased(int row) const
	{
		auto mask = WordConfiguration<Word>::AllOnes;
		if (row < WordConfiguration<Word>::WordSize-1) mask = ~(WordConfiguration<Word>::AllOnes << (row + 1));
		auto value = getScoreBeforeStart() + WordConfiguration<Word>::popcount(VP & mask) - WordConfiguration<Word>::popcount(VN & mask);
		return value;
	}

	ScoreType getPriorityScore(size_t j) const
	{
		ScoreType result = priorityScoreLocalMinima(j);
#ifdef EXTRACORRECTNESSASSERTIONS
		assert(result == priorityScoreCellByCell(j));
#endif
		return result;
	}

	double getChangedPriorityScore(WordSlice oldSlice, size_t j, double priorityMismatchPenalty) const
	{
		double result = changedPriorityScoreLocalMinima(oldSlice, j, priorityMismatchPenalty);
#ifdef EXTRACORRECTNESSASSERTIONS
		assert(result == changedPriorityScoreCellByCell(oldSlice, j, priorityMismatchPenalty));
#endif
		return result;
	}

	ScoreType getMinScore() const
	{
		ScoreType result = minScoreLocalMinima();
#ifdef EXTRACORRECTNESSASSERTIONS
		assert(result == minScoreCellByCell());
#endif
		return result;
	}

	ScoreType maxXScoreFirstSlices(ScoreType seqOffset, ScoreType errorCost, int cells, Word extraMask = WordConfiguration<Word>::AllOnes) const
	{
		assert(cells > 0);
		assert(cells <= WordConfiguration<Word>::WordSize);
		ScoreType result = maxXScoreLocalMinima(seqOffset, errorCost, cells, extraMask);
#ifdef EXTRACORRECTNESSASSERTIONS
		// todo figure out why breaks
		assert(result == maxXScoreCellByCell(seqOffset, errorCost, cells));
#endif
		return result;
	}

	ScoreType maxXScore(ScoreType seqOffset, ScoreType errorCost, Word extraMask = WordConfiguration<Word>::AllOnes) const
	{
		return maxXScoreFirstSlices(seqOffset, errorCost, WordConfiguration<Word>::WordSize, extraMask);
	}

	ScoreType getXScore(ScoreType seqOffset, int offset, ScoreType errorCost) const
	{
		return (seqOffset + offset+1)*100 - getValue(offset) * errorCost;
	}

	// ScoreType getScoreBeforeStart() const
	// {
	// 	return scoreEnd - WordConfiguration<Word>::popcount(VP) + WordConfiguration<Word>::popcount(VN);
	// }
	__attribute__((noinline))
	ScoreType getScoreBeforeStart() const
	{

		// --- Debug logging block ---
		if (enableGetScoreBeforeStartDebug) {
			getScoreBeforeStartIteration++;
			std::ofstream dbg("GbvCallTrace.log", std::ios::app);
			dbg << "getScoreBeforeStart call #" << getScoreBeforeStartIteration
				<< " | VP=" << VP
				<< " | VN=" << VN
				<< " | scoreEnd=" << scoreEnd
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
		// --- End debug logging block ---

		auto& regfile = GraphAlignerBitvectorCommon<LengthType, ScoreType, Word>::regfile;

		regfile[11] = scoreEnd;
		regfile[2] = VN;
		regfile[3] = VP;
		regfile[23] = WordConfiguration<Word>::popcount(VP);
		regfile[24] = WordConfiguration<Word>::popcount(VN);
		regfile[25] = regfile[11] - regfile[23];
		regfile[25] = regfile[25] + regfile[24];

		return regfile[25];
	}

#ifdef NDEBUG
	__attribute__((always_inline))
#endif
	ScoreType changedMinScore(WordSlice other) const
	{
		auto result = changedMinScoreLocalMinima(other);
#ifdef EXTRACORRECTNESSASSERTIONS
		assert(result == changedMinScoreCellByCell(other));
#endif
		return result;
	}

private:

#ifdef EXTRACORRECTNESSASSERTIONS
	ScoreType maxXScoreCellByCell(ScoreType seqOffset, ScoreType errorCost, size_t cells) const
	{
		ScoreType result = std::numeric_limits<ScoreType>::min();
		for (int i = 0; i < WordConfiguration<Word>::WordSize && i < cells; i++)
		{
			result = std::max(result, getXScore(seqOffset, i, errorCost));
		}
		return result;
	}
	ScoreType changedPriorityScoreCellByCell(WordSlice other, size_t j, double priorityMismatchPenalty) const
	{
		ScoreType result = std::numeric_limits<ScoreType>::max();
		if (getScoreBeforeStart() < other.getScoreBeforeStart()) result = getScoreBeforeStart()*priorityMismatchPenalty - (ScoreType)j + 1;
		for (size_t i = 0; i < WordConfiguration<Word>::WordSize; i++)
		{
			if (getValue(i) < other.getValue(i)) result = std::min(result, (ScoreType)(getValue(i)*priorityMismatchPenalty) - (ScoreType)j - (ScoreType)i);
		}
		return result;
	}
	ScoreType priorityScoreCellByCell(size_t j) const
	{
		ScoreType result = getScoreBeforeStart() - j / 2;
		for (size_t i = 0; i < WordConfiguration<Word>::WordSize; i++)
		{
			result = std::min(result, (ScoreType)(getValue(i) - (j+i+1)/2));
		}
		return result;
	}
	ScoreType changedMinScoreCellByCell(WordSlice other) const
	{
		ScoreType result = std::numeric_limits<ScoreType>::max();
		if (getScoreBeforeStart() < other.getScoreBeforeStart()) result = getScoreBeforeStart();
		for (size_t i = 0; i < WordConfiguration<Word>::WordSize; i++)
		{
			if (getValue(i) < other.getValue(i)) result = std::min(result, getValue(i));
		}
		return result;
	}
	ScoreType minScoreCellByCell() const
	{
		ScoreType minScore = std::numeric_limits<ScoreType>::max();
		for (int i = 0; i < WordConfiguration<Word>::WordSize; i++)
		{
			minScore = std::min(minScore, getValue(i));
		}
		return minScore;
	}
#endif

	ScoreType maxXScoreLocalMinima(ScoreType seqOffset, ScoreType errorCost, size_t cells, Word extraMask) const
	{
		ScoreType scoreBeforeStart = getScoreBeforeStart();
		//rightmost VP between any VN's, aka one cell to the left of a minimum
		Word priorityCausedMinima = ~VP;
		Word possibleLocalMinima = (VP & (priorityCausedMinima - VP));
		//shift right by one to get the minimum
		possibleLocalMinima >>= 1;
		//leftmost bit might be a minimum if there is no VP to its right
		possibleLocalMinima |= WordConfiguration<Word>::LastBit & (priorityCausedMinima | ~(priorityCausedMinima - VP)) & ~VP;
		ScoreType result = std::numeric_limits<ScoreType>::min();
		possibleLocalMinima &= extraMask;
		possibleLocalMinima |= 1;
		if (cells < WordConfiguration<Word>::WordSize)
		{
			possibleLocalMinima |= (Word)1 << (Word)(cells-1);
			possibleLocalMinima &= ~(WordConfiguration<Word>::AllOnes << cells);
		}
		else
		{
			possibleLocalMinima |= (Word)1 << (Word)(WordConfiguration<Word>::WordSize-1);
		}
		ScoreType zeroScore = seqOffset * 100 - scoreBeforeStart * errorCost;
		while (possibleLocalMinima != 0)
		{
			//all cells from the right up to the first minimum are one
			Word currentMinimumMask = possibleLocalMinima ^ (possibleLocalMinima-1);
			ScoreType cellsHere = WordConfiguration<Word>::popcount(currentMinimumMask);
			ScoreType scoreHere = (ScoreType)WordConfiguration<Word>::popcount(VP & currentMinimumMask) - (ScoreType)WordConfiguration<Word>::popcount(VN & currentMinimumMask);
			result = std::max(result, cellsHere*100 - scoreHere * errorCost);
			possibleLocalMinima &= ~currentMinimumMask;
		}
		result += zeroScore;
		return result;
	}

	double changedPriorityScoreLocalMinima(WordSlice oldSlice, size_t j, double priorityMismatchPenalty) const
	{
		ScoreType otherScoreBeforeStart = oldSlice.getScoreBeforeStart();
		ScoreType scoreBeforeStart = getScoreBeforeStart();
		//rightmost VP between any VN's, aka one cell to the left of a minimum
		Word priorityCausedMinima = ~VP;
		Word possibleLocalMinima = (VP & (priorityCausedMinima - VP));
		//shift right by one to get the minimum
		possibleLocalMinima >>= 1;
		//leftmost bit might be a minimum if there is no VP to its right
		possibleLocalMinima |= WordConfiguration<Word>::LastBit & (priorityCausedMinima | ~(priorityCausedMinima - VP)) & ~VP;
		if (scoreEnd + WordConfiguration<Word>::popcount(VN) >= oldSlice.scoreEnd - WordConfiguration<Word>::popcount(oldSlice.VP))
		{
			auto masks = differenceMasks(VP, VN, oldSlice.VP, oldSlice.VN, otherScoreBeforeStart - scoreBeforeStart);
			auto smaller = masks.first;
			//corner cases
			if (smaller != WordConfiguration<Word>::AllOnes)
			{
				possibleLocalMinima |= (~smaller) >> 1;
				possibleLocalMinima |= (~smaller) << 1;
				possibleLocalMinima |= 1;
				possibleLocalMinima |= WordConfiguration<Word>::LastBit;
				possibleLocalMinima &= smaller;
			}
		}
		ScoreType result = (scoreBeforeStart < otherScoreBeforeStart) ? (scoreBeforeStart*priorityMismatchPenalty-j+1) : std::numeric_limits<ScoreType>::max();
		while (possibleLocalMinima != 0)
		{
			//all cells from the right up to the first minimum are one
			Word currentMinimumMask = possibleLocalMinima ^ (possibleLocalMinima-1);
			ScoreType scoreHere = scoreBeforeStart + WordConfiguration<Word>::popcount(VP & currentMinimumMask) - WordConfiguration<Word>::popcount(VN & currentMinimumMask);
			scoreHere *= priorityMismatchPenalty;
			scoreHere -= (ScoreType)((j + WordConfiguration<Word>::popcount(currentMinimumMask))) - 1;
			result = std::min(result, scoreHere);
			possibleLocalMinima &= ~currentMinimumMask;
		}
		return result;
	}

	ScoreType priorityScoreLocalMinima(size_t j) const
	{
		ScoreType scoreBeforeStart = getScoreBeforeStart();
		//rightmost VP between any VN's, aka one cell to the left of a minimum
		Word priorityCausedMinima = 0xAAAAAAAAAAAAAAAA & ~VP & ~VN;
		priorityCausedMinima |= VN;
		Word possibleLocalMinima = (VP & (priorityCausedMinima - VP));
		//shift right by one to get the minimum
		possibleLocalMinima >>= 1;
		//leftmost bit might be a minimum if there is no VP to its right
		possibleLocalMinima |= WordConfiguration<Word>::LastBit & (priorityCausedMinima | ~(priorityCausedMinima - VP)) & ~VP;
		ScoreType result = scoreBeforeStart - j/2;
		while (possibleLocalMinima != 0)
		{
			//all cells from the right up to the first minimum are one
			Word currentMinimumMask = possibleLocalMinima ^ (possibleLocalMinima-1);
			ScoreType scoreHere = scoreBeforeStart + WordConfiguration<Word>::popcount(VP & currentMinimumMask) - WordConfiguration<Word>::popcount(VN & currentMinimumMask);
			scoreHere -= (ScoreType)((j + WordConfiguration<Word>::popcount(currentMinimumMask)) / 2);
			result = std::min(result, scoreHere);
			possibleLocalMinima &= ~currentMinimumMask;
		}
		return result;
	}

	ScoreType minScoreLocalMinima() const
	{
		ScoreType scoreBeforeStart = getScoreBeforeStart();
		//rightmost VP between any VN's, aka one cell to the left of a minimum
		Word possibleLocalMinima = (VP & (VN - VP));
		//shift right by one to get the minimum
		possibleLocalMinima >>= 1;
		//leftmost bit might be a minimum if there is no VP to its right
		possibleLocalMinima |= WordConfiguration<Word>::LastBit & (VN | ~(VN - VP)) & ~VP;
		ScoreType result = scoreBeforeStart + (VP & 1) - (VN & 1);
		//the score is inited to the first cell at the start
		possibleLocalMinima &= ~((Word)1);
		while (possibleLocalMinima != 0)
		{
			//all cells from the right up to the first minimum are one
			Word currentMinimumMask = possibleLocalMinima ^ (possibleLocalMinima-1);
			ScoreType scoreHere = scoreBeforeStart + WordConfiguration<Word>::popcount(VP & currentMinimumMask) - WordConfiguration<Word>::popcount(VN & currentMinimumMask);
			result = std::min(result, scoreHere);
			possibleLocalMinima &= ~currentMinimumMask;
		}
		return result;
	}

#ifdef NDEBUG
	__attribute__((always_inline))
#endif
	ScoreType changedMinScoreLocalMinima(WordSlice oldSlice) const
	{
		ScoreType scoreBeforeStart = getScoreBeforeStart();
		ScoreType otherScoreBeforeStart = oldSlice.getScoreBeforeStart();
		//rightmost VP between any VN's, aka one cell to the left of a minimum
		Word possibleLocalMinima = (VP & (VN - VP));
		//shift right by one to get the minimum
		possibleLocalMinima >>= 1;
		//leftmost bit might be a minimum if there is no VP to its right
		possibleLocalMinima |= WordConfiguration<Word>::LastBit & (VN | ~(VN - VP)) & ~VP;
		if (scoreEnd + WordConfiguration<Word>::popcount(VN) >= oldSlice.scoreEnd - WordConfiguration<Word>::popcount(oldSlice.VP))
		{
			auto masks = differenceMasks(VP, VN, oldSlice.VP, oldSlice.VN, otherScoreBeforeStart - scoreBeforeStart);
			auto smaller = masks.first;
			//corner cases
			if (smaller != WordConfiguration<Word>::AllOnes)
			{
				possibleLocalMinima |= (~smaller) >> 1;
				possibleLocalMinima |= (~smaller) << 1;
				possibleLocalMinima |= 1;
				possibleLocalMinima |= WordConfiguration<Word>::LastBit;
				possibleLocalMinima &= smaller;
			}
		}
		ScoreType result = (scoreBeforeStart < otherScoreBeforeStart) ? scoreBeforeStart : std::numeric_limits<ScoreType>::max();
		while (possibleLocalMinima != 0)
		{
			//all cells from the right up to the first minimum are one
			Word currentMinimumMask = possibleLocalMinima ^ (possibleLocalMinima-1);
			ScoreType scoreHere = scoreBeforeStart + WordConfiguration<Word>::popcount(VP & currentMinimumMask) - WordConfiguration<Word>::popcount(VN & currentMinimumMask);
			result = std::min(result, scoreHere);
			possibleLocalMinima &= ~currentMinimumMask;
		}
		return result;
	}

	static uint64_t bytePrefixSums(uint64_t value, int addition)
	{
		value <<= WordConfiguration<Word>::ChunkBits;
		assert(addition >= 0);
		value += addition;
		return value * WordConfiguration<Word>::PrefixSumMultiplierConstant;
	}

	static uint64_t bytePrefixSums(uint64_t value)
	{
		value <<= WordConfiguration<Word>::ChunkBits;
		return value * WordConfiguration<Word>::PrefixSumMultiplierConstant;
	}

	static uint64_t byteVPVNSum(uint64_t prefixSumVP, uint64_t prefixSumVN)
	{
		uint64_t result = WordConfiguration<Word>::SignMask;
		assert((prefixSumVP & result) == 0);
		assert((prefixSumVN & result) == 0);
		result += prefixSumVP;
		result -= prefixSumVN;
		result ^= WordConfiguration<Word>::SignMask;
		return result;
	}

// #ifndef NDEBUG
// 	__attribute__((always_inline))
// #endif
	// static WordSlice mergeTwoSlices(WordSlice left, WordSlice right)
	// {
	// 	//O(log w), because prefix sums need log w chunks of log w bits
	// 	static_assert(std::is_same<Word, uint64_t>::value);
	// 	if (left.getScoreBeforeStart() > right.getScoreBeforeStart()) std::swap(left, right);
	// 	assert((left.VP & left.VN) == WordConfiguration<Word>::AllZeros);
	// 	assert((right.VP & right.VN) == WordConfiguration<Word>::AllZeros);
	// 	auto masks = differenceMasks(left.VP, left.VN, right.VP, right.VN, right.getScoreBeforeStart() - left.getScoreBeforeStart());
	// 	return mergeTwoSlices(left, right, masks.first, masks.second);
	// }
	__attribute__((noinline))
	static WordSlice mergeTwoSlices(WordSlice left, WordSlice right)
	{

		// --- Debug logging block ---
		if (enableMergeTwoSlices2InputDebug) {
			mergeTwoSlices2InputIteration++;
			std::ofstream dbg("GbvCallTrace.log", std::ios::app);
			dbg << "mergeTwoSlices2Input call #" << mergeTwoSlices2InputIteration
				<< " | left.VP=" << left.VP
				<< " | left.VN=" << left.VN
				<< " | left.scoreEnd=" << left.scoreEnd
				<< " | right.VP=" << right.VP
				<< " | right.VN=" << right.VN
				<< " | right.scoreEnd=" << right.scoreEnd
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
		// --- End debug logging block ---

		auto& regfile = GraphAlignerBitvectorCommon<LengthType, ScoreType, Word>::regfile;
		static_assert(std::is_same<Word, uint64_t>::value);

		regfile[12] = left.VN;
		regfile[13] = left.VP;
		regfile[14] = left.getScoreBeforeStart();
		regfile[15] = left.scoreEnd;
		regfile[16] = right.VN;
		regfile[17] = right.VP;
		regfile[18] = right.getScoreBeforeStart();
		regfile[19] = right.scoreEnd;

		if (regfile[14] > regfile[18]) // std::swap(left, right);
		{ 
			regfile[22] = regfile[12];
			regfile[23] = regfile[13];
			regfile[24] = regfile[15];
			regfile[25] = regfile[14];

			regfile[12] = regfile[16];
			regfile[13] = regfile[17];
			regfile[15] = regfile[19];
			regfile[14] = regfile[18];
			
			regfile[16] = regfile[22];
			regfile[17] = regfile[23];
			regfile[19] = regfile[24];
			regfile[18] = regfile[25];
		}

		assert((regfile[13] & regfile[12]) == WordConfiguration<Word>::AllZeros);
		assert((regfile[17] & regfile[16]) == WordConfiguration<Word>::AllZeros);
		regfile[31] = regfile[18] - regfile[14];
		left.VN = regfile[12];
		left.VP = regfile[13];
		left.scoreEnd = regfile[15];
		right.VN = regfile[16];
		right.VP = regfile[17];
		right.scoreEnd = regfile[19];
		auto masks = differenceMasks(regfile[13], regfile[12], regfile[17], regfile[16], regfile[31]);
		regfile[20] = masks.first;
		regfile[21] = masks.second;

		return mergeTwoSlices(left, right, regfile[20], regfile[21]);
	}

// #ifdef NDEBUG
// 	__attribute__((always_inline))
// #endif
	// static WordSlice mergeTwoSlices(WordSlice left, WordSlice right, Word leftSmaller, Word rightSmaller)
	// {
	// 	assert(left.getScoreBeforeStart() <= right.getScoreBeforeStart());
	// 	WordSlice result;
	// 	assert((left.VP & left.VN) == WordConfiguration<Word>::AllZeros);
	// 	assert((right.VP & right.VN) == WordConfiguration<Word>::AllZeros);
	// 	assert((leftSmaller & rightSmaller) == 0);
	// 	auto mask = (rightSmaller | ((leftSmaller | rightSmaller) - (rightSmaller << 1))) & ~leftSmaller;
	// 	uint64_t leftReduction = leftSmaller & (rightSmaller << 1);
	// 	uint64_t rightReduction = rightSmaller & (leftSmaller << 1);
	// 	if ((rightSmaller & 1) && left.getScoreBeforeStart() < right.getScoreBeforeStart())
	// 	{
	// 		rightReduction |= 1;
	// 	}
	// 	assert((leftReduction & right.VP) == leftReduction);
	// 	assert((rightReduction & left.VP) == rightReduction);
	// 	assert((leftReduction & left.VN) == leftReduction);
	// 	assert((rightReduction & right.VN) == rightReduction);
	// 	left.VN &= ~leftReduction;
	// 	right.VN &= ~rightReduction;
	// 	result.VN = (left.VN & ~mask) | (right.VN & mask);
	// 	result.VP = (left.VP & ~mask) | (right.VP & mask);
	// 	assert((result.VP & result.VN) == 0);
	// 	result.scoreEnd = std::min(left.scoreEnd, right.scoreEnd);
	// 	return result;
	// }
	
	__attribute__((noinline))
	static WordSlice mergeTwoSlices(WordSlice left, WordSlice right, Word leftSmaller, Word rightSmaller)
	{

		// --- Debug logging block ---
		if (enableMergeTwoSlices4InputDebug) {
			mergeTwoSlices4InputIteration++;
			std::ofstream dbg("GbvCallTrace.log", std::ios::app);
			dbg << "mergeTwoSlices4Input call #" << mergeTwoSlices4InputIteration
				<< " | left.VP=" << left.VP
				<< " | left.VN=" << left.VN
				<< " | left.scoreEnd=" << left.scoreEnd
				<< " | right.VP=" << right.VP
				<< " | right.VN=" << right.VN
				<< " | right.scoreEnd=" << right.scoreEnd
				<< " | leftSmaller=" << leftSmaller
				<< " | rightSmaller=" << rightSmaller
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
		// --- End debug logging block ---

		auto& regfile = GraphAlignerBitvectorCommon<LengthType, ScoreType, Word>::regfile;

		WordSlice result;

		regfile[12] = left.VN;
		regfile[13] = left.VP;
		regfile[14] = left.getScoreBeforeStart();
		regfile[15] = left.scoreEnd;
		regfile[16] = right.VN;
		regfile[17] = right.VP;
		regfile[18] = right.getScoreBeforeStart();
		regfile[19] = right.scoreEnd;

		regfile[28] = result.VN;
		regfile[29] = result.VP;
		regfile[30] = result.scoreEnd;

		regfile[20] = leftSmaller;
		regfile[21] = rightSmaller;

		assert(regfile[14] <= regfile[18]);
		assert((regfile[13] & regfile[12]) == WordConfiguration<Word>::AllZeros);

		assert((regfile[17] & regfile[16]) == WordConfiguration<Word>::AllZeros);

/* 		static int assert_print_count = 0;
		if (assert_print_count < 5) {
            std::cerr << "Assertion Check " << assert_print_count + 1 << "\n";
			std::cerr << "regfile[20] = 0x" << std::hex << regfile[20] << " (" << std::dec << regfile[20] << ")\n";
			std::cerr << "regfile[21] = 0x" << std::hex << regfile[21] << " (" << std::dec << regfile[21] << ")\n";
			std::cerr << "regfile[20] & regfile[21] = 0x" << std::hex << (regfile[20] & regfile[21]) << " (" << std::dec << (regfile[20] & regfile[21]) << ")\n";
			assert_print_count++;
} */
		assert((regfile[20] & regfile[21]) == 0);

		regfile[22] = regfile[20] | regfile[21]; // leftsmaller | rightsmaller
		regfile[23] = regfile[21] << 1; // rightsmaller << 1
		regfile[22] = regfile[22] - regfile[23]; // ((leftSmaller | rightSmaller) - (rightSmaller << 1))
		regfile[22] = regfile[21] | regfile[22]; // rightsmaller | regfile25
		regfile[23] = ~regfile[20]; // ~leftsmaller
		regfile[22] = regfile[22] & regfile[23]; // reg23 = mask

		regfile[23] = regfile[21] << 1; // rightsmaller << 1
		regfile[23] = regfile[20] & regfile[23]; // leftreduction = reg24
		regfile[24] = regfile[20] << 1; // leftsmaller << 1
		regfile[24] = regfile[21] & regfile[24]; // rightreduction = reg25

		// if ((regfile[21] & 1) && regfile[14] < regfile[18])
		// {
		// 	regfile[24] |= 1;
		// }

		// regfile[24] = ((regfile[21] & 1) && (regfile[14] < regfile[18])) ? (regfile[24] | 1) : regfile[24];
		regfile[0] = 0;
		regfile[25] = 1;
		regfile[26] = regfile[21] & regfile[25];
		regfile[27] = (regfile[18] > regfile[14]) ? regfile[25] : regfile[0];
		regfile[28] = regfile[26] & regfile[27];
		regfile[26] = regfile[24] | regfile[25];
		regfile[24] = (regfile[28] == regfile[25]) ? regfile[26] : regfile[24];

		assert((regfile[23] & regfile[17]) == regfile[23]); 
		assert((regfile[24] & regfile[13]) == regfile[24]);  
		assert((regfile[23] & regfile[12]) == regfile[23]);  
		assert((regfile[24] & regfile[16]) == regfile[24]); 

		regfile[25] = ~regfile[23]; // ~leftreduction
		regfile[12] = regfile[12] & regfile[25]; // leftvn &= ~leftreduction
		regfile[25] = ~regfile[24]; // ~rightreduction
		regfile[16] = regfile[16] & regfile[25]; // rightvn &= ~rightreduction

		regfile[25] = ~regfile[22]; //~mask
		regfile[26] = regfile[12] & regfile[25]; // leftVN & ~mask
		regfile[27] = regfile[16] & regfile[22]; // rightVN & mask
		regfile[28] = regfile[26] | regfile[27]; // (left.VN & ~mask) | (right.VN & mask);

		regfile[26] = regfile[13] & regfile[25]; // leftVP & ~mask
		regfile[27] = regfile[17] & regfile[22]; // rightVP & mask
		regfile[29] = regfile[26] | regfile[27]; // (left.VP & ~mask) | (right.VP & mask);

		assert((regfile[28] & regfile[29]) == 0);

		regfile[30] = std::min(regfile[15], regfile[19]);
		
		result.VN = regfile[28];
		result.VP = regfile[29];
		result.scoreEnd = regfile[30];
		left.VN = regfile[12];
		right.VN = regfile[16];

		return result;
	}

#ifdef NDEBUG
	__attribute__((always_inline))
#endif
	static std::pair<uint64_t, uint64_t> differenceMasks(uint64_t leftVP, uint64_t leftVN, uint64_t rightVP, uint64_t rightVN, int scoreDifference)
	{
		
		auto result = differenceMasksBitTwiddle(leftVP, leftVN, rightVP, rightVN, scoreDifference);
#ifdef EXTRACORRECTNESSASSERTIONS
		auto debugCompare = differenceMasksWord(leftVP, leftVN, rightVP, rightVN, scoreDifference);
		assert(result.first == debugCompare.first);
		assert(result.second == debugCompare.second);
#endif
		return result;
	}

	static ScoreType clamp(ScoreType low, ScoreType val, ScoreType high)
	{
		return std::min(high, std::max(low, val));
	}


/*
#ifdef NDEBUG
	__attribute__((always_inline))
#endif
	__attribute__((optimize("unroll-loops")))
	static std::pair<Word, Word> differenceMasksBitTwiddle(Word leftVP, Word leftVN, Word rightVP, Word rightVN, int scoreDifference)
	{
		Word leftSmaller = 0;
		Word rightSmaller = 0;

			static int print_count_5 = 0;
	if (print_count_5 < 5) {
		std::cerr << "initial values " << print_count_5 + 1 << "\n";
		std::cerr << "  leftVP: 0x" << std::hex << leftVP << " (" << std::dec << leftVP << ")\n";
		std::cerr << "  rightVP: 0x" << std::hex << rightVP << " (" << std::dec << rightVP << ")\n";
		std::cerr << "  leftVN: 0x" << std::hex << leftVN << " (" << std::dec << leftVN << ")\n";
		std::cerr << "  rightVN: 0x" << std::hex << rightVN << " (" << std::dec << rightVN << ")\n";
		print_count_5++;
	}

		Word VPcommon = ~(leftVP & rightVP);
		Word VNcommon = ~(leftVN & rightVN);
		leftVP &= VPcommon;
		leftVN &= VNcommon;
		rightVP &= VPcommon;
		rightVN &= VNcommon;

	static int print_count = 0;
	if (print_count < 5) {
		std::cerr << "Basic values " << print_count + 1 << "\n";
		std::cerr << "  VPCommon: 0x" << std::hex << VPcommon << " (" << std::dec << VPcommon << ")\n";
		std::cerr << "  VNCommon: 0x" << std::hex << VNcommon << " (" << std::dec << VNcommon << ")\n";
		std::cerr << "  leftVP: 0x" << std::hex << leftVP << " (" << std::dec << leftVP << ")\n";
		std::cerr << "  rightVP: 0x" << std::hex << rightVP << " (" << std::dec << rightVP << ")\n";
		std::cerr << "  leftVN: 0x" << std::hex << leftVN << " (" << std::dec << leftVN << ")\n";
		std::cerr << "  rightVN: 0x" << std::hex << rightVN << " (" << std::dec << rightVN << ")\n";
		print_count++;
	}

		Word twosmaller = leftVN & rightVP; //left is two smaller
		Word onesmaller = (rightVP & ~leftVN) | (leftVN & ~rightVP);
		Word onebigger = (leftVP & ~rightVN) | (rightVN & ~leftVP);
		Word twobigger = rightVN & leftVP; //left is two bigger
		onebigger |= twobigger;
		onesmaller |= twosmaller;

					// DEBUG
			static int print_count_1 = 0;
			if (print_count_1 < 5) {
            std::cerr << "diff Mask Initial " << print_count_1 + 1 << "\n";
            std::cerr << "  leftSmaller: 0x" << std::hex << leftSmaller << " (" << std::dec << leftSmaller << ")\n";
			std::cerr << "  rightSmaller: 0x" << std::hex << rightSmaller << " (" << std::dec << rightSmaller << ")\n";
            std::cerr << "  onebigger: 0x" << std::hex << onebigger << " (" << std::dec << onebigger << ")\n";
            std::cerr << "  twobigger: 0x" << std::hex << twobigger << " (" << std::dec << twobigger << ")\n";
            std::cerr << "  onesmaller: 0x" << std::hex << onesmaller << " (" << std::dec << onesmaller << ")\n";
            std::cerr << "  twosmaller: 0x" << std::hex << twosmaller << " (" << std::dec << twosmaller << ")\n";
            print_count_1++;
        	}
			//END DEBUG

		//scoredifference is right - left
		if (scoreDifference > 0)
		{
			//right is higher
			for (int i = 1; i < scoreDifference; i++)
			{
				Word leastSignificant = onebigger & ~(onebigger - 1);
				onebigger ^= (~twobigger & leastSignificant);
				twobigger &= ~leastSignificant;

				if (onebigger == 0)
				{
					return std::make_pair(WordConfiguration<Word>::AllOnes, WordConfiguration<Word>::AllZeros);
				}
			}

			Word leastSignificant = onebigger & ~(onebigger - 1);
			leftSmaller |= leastSignificant - 1;
			onebigger ^= (~twobigger & leastSignificant);
			twobigger &= ~leastSignificant;


		}
		else if (scoreDifference < 0)
		{
			//left is higher
			for (int i = 1; i < -scoreDifference; i++)
			{
				Word leastSignificant = onesmaller & ~(onesmaller - 1);
				onesmaller ^= (~twosmaller & leastSignificant);
				twosmaller &= ~leastSignificant;
				if (onesmaller == 0)
				{
					return std::make_pair(WordConfiguration<Word>::AllZeros, WordConfiguration<Word>::AllOnes);
				}
			}
			Word leastSignificant = onesmaller & ~(onesmaller - 1);
			rightSmaller |= leastSignificant - 1;
			onesmaller ^= (~twosmaller & leastSignificant);
			twosmaller &= ~leastSignificant;

						// DEBUG
			static int print_count_2 = 0;
			if (print_count_2 < 5) {
            std::cerr << "[differenceMasksBitTwiddle] Iteration " << print_count_2 + 1 << "\n";
            std::cerr << "  rightSmaller: 0x" << std::hex << rightSmaller << " (" << std::dec << rightSmaller << ")\n";
            std::cerr << "  leastSignificant: 0x" << std::hex << leastSignificant << " (" << std::dec << leastSignificant << ")\n";
            std::cerr << "  onebigger: 0x" << std::hex << onebigger << " (" << std::dec << onebigger << ")\n";
            std::cerr << "  twobigger: 0x" << std::hex << twobigger << " (" << std::dec << twobigger << ")\n";
            std::cerr << "  onesmaller: 0x" << std::hex << onesmaller << " (" << std::dec << onesmaller << ")\n";
            std::cerr << "  twosmaller: 0x" << std::hex << twosmaller << " (" << std::dec << twosmaller << ")\n";
            print_count_2++;
        	}
			//END DEBUG
		}
		for (int i = 0; i < WordConfiguration<Word>::WordSize; i++)
		{
			if (onesmaller == 0)
			{
				if (onebigger == 0) break;
											// DEBUG
			static int print_count_11 = 0;
			if (print_count_11 < 5) {
            std::cerr << " Line 870 initial " << print_count_11 + 1 << "\n";
            std::cerr << "  rightSmaller: 0x" << std::hex << rightSmaller << " (" << std::dec << rightSmaller << ")\n";
            print_count_11++;
        	}
			//END DEBUG

				Word leastSignificant = onebigger & ~(onebigger - 1);
				rightSmaller |= -leastSignificant;

							// DEBUG
			static int print_count_3 = 0;
			if (print_count_3 < 5) {
            std::cerr << " Line 870 " << print_count_3 + 1 << "\n";
            std::cerr << "  rightSmaller: 0x" << std::hex << rightSmaller << " (" << std::dec << rightSmaller << ")\n";
            std::cerr << "  leastSignificant: 0x" << std::hex << leastSignificant << " (" << std::dec << leastSignificant << ")\n";
            std::cerr << "  -leastSignificant: 0x" << std::hex << -leastSignificant << " (" << std::dec << -leastSignificant << ")\n";
            std::cerr << "  onebigger: 0x" << std::hex << onebigger << " (" << std::dec << onebigger << ")\n";
            std::cerr << "  twobigger: 0x" << std::hex << twobigger << " (" << std::dec << twobigger << ")\n";
            std::cerr << "  onesmaller: 0x" << std::hex << onesmaller << " (" << std::dec << onesmaller << ")\n";
            std::cerr << "  twosmaller: 0x" << std::hex << twosmaller << " (" << std::dec << twosmaller << ")\n";
            print_count_3++;
        	}
			//END DEBUG
				break;
			}
			if (onebigger == 0)
			{
#ifdef EXTRACORRECTNESSASSERTIONS
				assert(onesmaller != 0);
#endif
				Word leastSignificant = onesmaller & ~(onesmaller - 1);
				leftSmaller |= -leastSignificant;


				break;
			}
			Word leastSignificantBigger = onebigger & ~(onebigger - 1);
			Word leastSignificantSmaller = onesmaller & ~(onesmaller - 1);
#ifdef EXTRACORRECTNESSASSERTIONS
			assert((onebigger & leastSignificantBigger) != 0);
			assert((onesmaller & leastSignificantSmaller) != 0);
			assert(leastSignificantSmaller != leastSignificantBigger);
			assert(leastSignificantSmaller != 0);
			assert(leastSignificantBigger != 0);
#endif
			if (leastSignificantBigger > leastSignificantSmaller)
			{
				leftSmaller |= leastSignificantBigger - leastSignificantSmaller;

			}
			else
			{
				rightSmaller |= leastSignificantSmaller - leastSignificantBigger;

				
							// DEBUG
			static int print_count_4 = 0;
			if (print_count_4 < 100000) {
            std::cerr << "[differenceMasksBitTwiddle] Line 900 " << print_count_4 + 1 << "\n";
            std::cerr << "  rightSmaller: 0x" << std::hex << rightSmaller << " (" << std::dec << rightSmaller << ")\n";
			std::cerr << "  leastSignificantBigger: 0x" << std::hex << leastSignificantBigger << " (" << std::dec << leastSignificantBigger << ")\n";
			std::cerr << "  leastSignificantSmaller: 0x" << std::hex << leastSignificantSmaller << " (" << std::dec << leastSignificantSmaller << ")\n";
            std::cerr << "  onebigger: 0x" << std::hex << onebigger << " (" << std::dec << onebigger << ")\n";
            std::cerr << "  twobigger: 0x" << std::hex << twobigger << " (" << std::dec << twobigger << ")\n";
            std::cerr << "  onesmaller: 0x" << std::hex << onesmaller << " (" << std::dec << onesmaller << ")\n";
            std::cerr << "  twosmaller: 0x" << std::hex << twosmaller << " (" << std::dec << twosmaller << ")\n";
            print_count_4++;
        	}
			//END DEBUG
			}
			onebigger ^= (~twobigger & leastSignificantBigger);
			twobigger &= ~leastSignificantBigger;
			onesmaller ^= (~twosmaller & leastSignificantSmaller);
			twosmaller &= ~leastSignificantSmaller;
		}
#ifdef EXTRACORRECTNESSASSERTIONS
		assert((leftSmaller & rightSmaller) == 0);
		assert(onesmaller == 0 || onebigger == 0);
#endif
		return std::make_pair(leftSmaller, rightSmaller);
	} 
*/


#ifdef NDEBUG
	__attribute__((always_inline))
#endif
	__attribute__((optimize("unroll-loops")))
	static std::pair<Word, Word> differenceMasksBitTwiddle(Word leftVP, Word leftVN, Word rightVP, Word rightVN, int scoreDifference)
	{

		// --- Debug logging block ---
		if (enableDifferenceMasksBitTwiddleDebug) {
			differenceMasksBitTwiddleIteration++;
			std::ofstream dbg("GbvCallTrace.log", std::ios::app);
			dbg << "differenceMasks call #" << differenceMasksBitTwiddleIteration
				<< " | leftVP=" << leftVP
				<< " | leftVN=" << leftVN
				<< " | rightVP=" << rightVP
				<< " | rightVN=" << rightVN
				<< " | scoreDifference=" << scoreDifference
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
		// --- End debug logging block ---

		auto& regfile = GraphAlignerBitvectorCommon<LengthType, ScoreType, Word>::regfile;
		regfile[12] = leftVN;
		regfile[13] = leftVP;
		regfile[16] = rightVN;
		regfile[17] = rightVP;
		regfile[22] = static_cast<Word>(scoreDifference); // DVS: This is an integer, may cause problems for later. Be careful. 

		regfile[20] = 0; // leftSmaller
		regfile[21] = 0; // rightSmaller

		regfile[23] = regfile[13] & regfile[17]; // VPcommon
		regfile[23] = ~regfile[23];

		regfile[24] = regfile[12] & regfile[16]; // VNcommon
		regfile[24] = ~regfile[24];

		regfile[13] = regfile[13] & regfile[23]; // leftVP &= VPcommon;
		regfile[12] = regfile[12] & regfile[24]; // leftVN &= VNcommon;

		regfile[17] = regfile[17] & regfile[23]; // rightVP &= VPcommon;
		regfile[16] = regfile[16] & regfile[24]; // rightVN &= VNcommon;

		regfile[25] = regfile[12] & regfile[17]; // Word twosmaller = leftVN & rightVP; //left is two smaller

		// Word onesmaller = (rightVP & ~leftVN) | (leftVN & ~rightVP);
		regfile[26] = ~regfile[12]; // ~leftVN
		regfile[26] = regfile[17] & regfile[26]; // rightVP & ~leftVN
		regfile[27] = ~regfile[17]; // ~rightVP
		regfile[27] = regfile[12] & regfile[27]; // leftVN & ~rightVP
		regfile[26] = regfile[26] | regfile[27]; // onesmaller

		// Word onebigger = (leftVP & ~rightVN) | (rightVN & ~leftVP);
		regfile[27] = ~regfile[16];
		regfile[27] = regfile[13] & regfile[27];
		regfile[28] = ~regfile[13];
		regfile[28] = regfile[28] & regfile[16];
		regfile[27] = regfile[27] | regfile[28]; // onebigger

		regfile[28] = regfile[16] & regfile[13]; // Word twobigger = rightVN & leftVP; //left is two bigger
		
		regfile[27] = regfile[27] | regfile[28]; // onebigger |= twobigger;
		regfile[26] = regfile[26] | regfile[25]; // onesmaller |= twosmaller;
		
		//scoredifference is right - left
		if (regfile[22] > 0)
		{
			//right is higher
			for (int i = 1; i < regfile[22]; i++)
			{
				regfile[23] = regfile[27] - 1;
				regfile[23] = ~regfile[23];
				regfile[23] = regfile[23] & regfile[27]; // Word leastSignificant = onebigger & ~(onebigger - 1);
				
				regfile[24] = ~regfile[28];
				regfile[24] = regfile[24] & regfile[23];
				regfile[27] = regfile[27] ^ regfile[24]; // onebigger ^= (~twobigger & leastSignificant);

				regfile[24] = ~regfile[23];
				regfile[28] = regfile[28] & regfile[24]; // twobigger &= ~leastSignificant;
								
				if (regfile[27] == 0)
				{
					return std::make_pair(WordConfiguration<Word>::AllOnes, WordConfiguration<Word>::AllZeros);
				}
			}

			regfile[23] = regfile[27] - 1;
			regfile[23] = ~regfile[23];
			regfile[23] = regfile[23] & regfile[27]; // Word leastSignificant = onebigger & ~(onebigger - 1);
			
			regfile[24] = regfile[23] - 1;
			regfile[20] = regfile[20] | regfile[24]; // leftSmaller |= leastSignificant - 1;

			regfile[24] = ~regfile[28];
			regfile[24] = regfile[24] & regfile[23];
			regfile[27] = regfile[27] ^ regfile[24]; // onebigger ^= (~twobigger & leastSignificant);

			regfile[24] = ~regfile[23];
			regfile[28] = regfile[28] & regfile[24]; // twobigger &= ~leastSignificant;

		}

		else if (regfile[22] < 0)
		{
			//left is higher
			for (int i = 1; i < -scoreDifference; i++)
			{
				regfile[23] = regfile[26] - 1;
				regfile[23] = ~regfile[23];
				regfile[23] = regfile[23] & regfile[26]; // Word leastSignificant = onesmaller & ~(onesmaller - 1);

				regfile[24] = ~regfile[25];
				regfile[24] = regfile[24] & regfile[23];
				regfile[26] = regfile[26] ^ regfile[24]; // onesmaller ^= (~twosmaller & leastSignificant);

				regfile[24] = ~regfile[23];
				regfile[25] = regfile[25] & regfile[24]; // twosmaller &= ~leastSignificant;				
				
				if (regfile[26] == 0)
				{
					return std::make_pair(WordConfiguration<Word>::AllZeros, WordConfiguration<Word>::AllOnes);
				}
			}
			regfile[23] = regfile[26] - 1;
			regfile[23] = ~regfile[23];
			regfile[23] = regfile[23] & regfile[26]; // Word leastSignificant = onesmaller & ~(onesmaller - 1);
			
			regfile[24] = regfile[23] - 1;
			regfile[21] = regfile[24] | regfile[21]; // rightSmaller |= leastSignificant - 1;

			regfile[24] = ~regfile[25];
			regfile[24] = regfile[24] & regfile[23];
			regfile[26] = regfile[26] ^ regfile[24]; // onesmaller ^= (~twosmaller & leastSignificant);

			regfile[24] = ~regfile[23];
			regfile[28] = regfile[28] & regfile[24]; // twobigger &= ~leastSignificant;

/* 						// DEBUG
			static int print_count_2 = 0;
			if (print_count_2 < 5) {
            std::cerr << "[differenceMasksBitTwiddle] Iteration " << print_count_2 + 1 << "\n";
            std::cerr << "  rightSmaller: 0x" << std::hex << regfile[21] << " (" << std::dec << regfile[21] << ")\n";
            std::cerr << "  leastSignificant: 0x" << std::hex << regfile[23] << " (" << std::dec << regfile[23] << ")\n";
            std::cerr << "  onebigger: 0x" << std::hex << regfile[27] << " (" << std::dec << regfile[27] << ")\n";
            std::cerr << "  twobigger: 0x" << std::hex << regfile[28] << " (" << std::dec << regfile[28] << ")\n";
            std::cerr << "  onesmaller: 0x" << std::hex << regfile[26] << " (" << std::dec << regfile[26] << ")\n";
            std::cerr << "  twosmaller: 0x" << std::hex << regfile[25] << " (" << std::dec << regfile[25] << ")\n";
            print_count_2++;
        	}
			//END DEBUG */
		}
		
		for (int i = 0; i < WordConfiguration<Word>::WordSize; i++)
		{
			if (regfile[26] == 0)
			{
				if (regfile[27] == 0) break;
				regfile[23] = regfile[27] - 1;
				regfile[23] = ~regfile[23];
				regfile[23] = regfile[23] & regfile[27]; // Word leastSignificant = onebigger & ~(onebigger - 1);
				
				// regfile[24] = regfile[23]; // removed the negation
				regfile[21] = regfile[21] | regfile[23]; // rightSmaller |= -leastSignificant;
				break;
			}
			if (regfile[27] == 0)
			{
#ifdef EXTRACORRECTNESSASSERTIONS
				assert(onesmaller != 0);
#endif
				regfile[23] = regfile[26] - 1;
				regfile[23] = ~regfile[23];
				regfile[23] = regfile[23] & regfile[26]; // Word leastSignificant = onesmaller & ~(onesmaller - 1);

				// regfile[24] = regfile[23]; // Removed the negation
				regfile[20] = regfile[20] | regfile[23]; // leftSmaller |= -leastSignificant;
				break;
			}
			regfile[29] = regfile[27] - 1;
			regfile[29] = ~regfile[29];
			regfile[29] = regfile[29] & regfile[27]; // Word leastSignificantBigger = onebigger & ~(onebigger - 1);

			regfile[30] = regfile[26] - 1;
			regfile[30] = ~regfile[30];
			regfile[30] = regfile[30] & regfile[26];// Word leastSignificantSmaller = onesmaller & ~(onesmaller - 1);

#ifdef EXTRACORRECTNESSASSERTIONS
			assert((onebigger & leastSignificantBigger) != 0);
			assert((onesmaller & leastSignificantSmaller) != 0);
			assert(leastSignificantSmaller != leastSignificantBigger);
			assert(leastSignificantSmaller != 0);
			assert(leastSignificantBigger != 0);
#endif
			if (regfile[29] > regfile[30])
			{
				regfile[24] = regfile[29] - regfile[30];
				regfile[20] = regfile[20] | regfile[24]; // leftSmaller |= leastSignificantBigger - leastSignificantSmaller;
			}
			else
			{
				regfile[24] = regfile[30] - regfile[29];
				regfile[21] = regfile[21] | regfile[24]; // rightSmaller |= leastSignificantSmaller - leastSignificantBigger;


			}
			regfile[24] = ~regfile[28];
			regfile[24] = regfile[29] & regfile[24];
			regfile[27] = regfile[27] ^ regfile[24]; // onebigger ^= (~twobigger & leastSignificantBigger);

			regfile[24] = ~regfile[29];
			regfile[28] = regfile[28] & regfile[24]; //	twobigger &= ~leastSignificantBigger;

			regfile[24] = ~regfile[25];
			regfile[24] = regfile[24] & regfile[30];
			regfile[26] = regfile[24] ^ regfile[26]; //	onesmaller ^= (~twosmaller & leastSignificantSmaller);

			regfile[24] = ~regfile[30];
			regfile[25] = regfile[25] & regfile[24]; //	twosmaller &= ~leastSignificantSmaller;
		}
#ifdef EXTRACORRECTNESSASSERTIONS
		assert((leftSmaller & rightSmaller) == 0);
		assert(onesmaller == 0 || onebigger == 0);
#endif

		return std::make_pair(regfile[20], regfile[21]);
	} 
	



	static std::pair<uint64_t, uint64_t> differenceMasksWord(uint64_t leftVP, uint64_t leftVN, uint64_t rightVP, uint64_t rightVN, int scoreDifference)
	{
		assert(scoreDifference >= 0);
		const uint64_t signmask = WordConfiguration<Word>::SignMask;
		const uint64_t lsbmask = WordConfiguration<Word>::LSBMask;
		const int chunksize = WordConfiguration<Word>::ChunkBits;
		const uint64_t allones = WordConfiguration<Word>::AllOnes;
		const uint64_t allzeros = WordConfiguration<Word>::AllZeros;
		uint64_t VPcommon = ~(leftVP & rightVP);
		uint64_t VNcommon = ~(leftVN & rightVN);
		leftVP &= VPcommon;
		leftVN &= VNcommon;
		rightVP &= VPcommon;
		rightVN &= VNcommon;
		//left is lower everywhere
		if (scoreDifference > WordConfiguration<Word>::popcount(rightVN) + WordConfiguration<Word>::popcount(leftVP))
		{
			return std::make_pair(allones, allzeros);
		}
		if (scoreDifference == 128 && rightVN == allones && leftVP == allones)
		{
			return std::make_pair(allones ^ ((Word)1 << (WordConfiguration<Word>::WordSize-1)), allzeros);
		}
		else if (scoreDifference == 0 && rightVN == allones && leftVP == allones)
		{
			return std::make_pair(0, allones);
		}
		assert(scoreDifference >= 0);
		assert(scoreDifference < 128);
		uint64_t byteVPVNSumLeft = byteVPVNSum(bytePrefixSums(WordConfiguration<Word>::ChunkPopcounts(leftVP), 0), bytePrefixSums(WordConfiguration<Word>::ChunkPopcounts(leftVN), 0));
		uint64_t byteVPVNSumRight = byteVPVNSum(bytePrefixSums(WordConfiguration<Word>::ChunkPopcounts(rightVP), scoreDifference), bytePrefixSums(WordConfiguration<Word>::ChunkPopcounts(rightVN), 0));
		uint64_t difference = byteVPVNSumLeft;
		{
			//take the bytvpvnsumright and split it from positive/negative values into two vectors with positive values, one which needs to be added and the other deducted
			//smearmask is 1 where the number needs to be deducted, and 0 where it needs to be added
			//except sign bits which are all 0
			uint64_t smearmask = ((byteVPVNSumRight & signmask) >> (chunksize-1)) * ((((Word)1) << (chunksize-1))-1);
			assert((smearmask & signmask) == 0);
			uint64_t deductions = ~smearmask & byteVPVNSumRight & ~signmask;
			//byteVPVNSumRight is in one's complement so take the not-value + 1
			uint64_t additions = (smearmask & ~byteVPVNSumRight) + (smearmask & lsbmask);
			assert((deductions & signmask) == 0);
			uint64_t signsBefore = difference & signmask;
			//unset the sign bits so additions don't interfere with other chunks
			difference &= ~signmask;
			difference += additions;
			//the sign bit is 1 if the value went from <0 to >=0
			//so in that case we need to flip it
			difference ^= signsBefore;
			signsBefore = difference & signmask;
			//set the sign bits so that deductions don't interfere with other chunks
			difference |= signmask;
			difference -= deductions;
			//sign bit is 0 if the value went from >=0 to <0
			//so flip them to the correct values
			signsBefore ^= signmask & ~difference;
			difference &= ~signmask;
			difference |= signsBefore;
		}
		//difference now contains the prefix sum difference (left-right) at each chunk
		uint64_t resultLeftSmallerThanRight = 0;
		uint64_t resultRightSmallerThanLeft = 0;
		for (int bit = 0; bit < chunksize; bit++)
		{
			uint64_t signsBefore = difference & signmask;
			//unset the sign bits so additions don't interfere with other chunks
			difference &= ~signmask;
			difference += leftVP & lsbmask;
			difference += rightVN & lsbmask;
			//the sign bit is 1 if the value went from <0 to >=0
			//so in that case we need to flip it
			difference ^= signsBefore;
			signsBefore = difference & signmask;
			//set the sign bits so that deductions don't interfere with other chunks
			difference |= signmask;
			difference -= leftVN & lsbmask;
			difference -= rightVP & lsbmask;
			//sign bit is 0 if the value went from >=0 to <0
			//so flip them to the correct values
			signsBefore ^= signmask & ~difference;
			difference &= ~signmask;
			difference |= signsBefore;
			leftVN >>= 1;
			leftVP >>= 1;
			rightVN >>= 1;
			rightVP >>= 1;
			//difference now contains the prefix sums difference (left-right) at each byte at (bit)'th bit
			//left < right when the prefix sum difference is negative (sign bit is set)
			uint64_t negative = (difference & signmask);
			resultLeftSmallerThanRight |= negative >> (WordConfiguration<Word>::ChunkBits - 1 - bit);
			//Test equality to zero. If it's zero, substracting one will make the sign bit 0, otherwise 1
			uint64_t notEqualToZero = ((difference | signmask) - lsbmask) & signmask;
			//right > left when the prefix sum difference is positive (not zero and not negative)
			resultRightSmallerThanLeft |= (notEqualToZero & ~negative) >> (WordConfiguration<Word>::ChunkBits - 1 - bit);
		}
		return std::make_pair(resultLeftSmallerThanRight, resultRightSmallerThanLeft);
	}

};

#endif