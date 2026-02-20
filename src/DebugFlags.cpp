#include "DebugFlags.h"
#include <fstream>
#include <sstream>

bool enableCalculateSliceDebug = true;
uint64_t calculateSliceIteration = 0;

bool enableGetNextSliceDebug = true;
uint64_t getNextSliceIteration = 0;

bool enableCalculateNodeClipPreciseDebug = true;
uint64_t calculateNodeClipPreciseIteration = 0;

bool enableCalculateNodeInnerDebug = true;
uint64_t calculateNodeInnerIteration = 0;

bool enableGetScoreBeforeStartDebug = true;
uint64_t getScoreBeforeStartIteration = 0;

bool enableMergeTwoSlices2InputDebug = true;
uint64_t mergeTwoSlices2InputIteration = 0;

bool enableMergeTwoSlices4InputDebug = true;
uint64_t mergeTwoSlices4InputIteration = 0;

bool enableDifferenceMasksBitTwiddleDebug = true;
uint64_t differenceMasksBitTwiddleIteration = 0;

bool enableFlattenWordSliceDebug = true;
uint64_t flattenWordSliceIteration = 0;

uint64_t EqVectorIteration = 0;


bool debugTop = true;

// Helper function for debug logging
void DEBUG_LOG(const std::string& functionName,
               bool& enableFlag,
               uint64_t& iterationCounter,
               std::function<void(std::ostream&)> outputFunc,
               bool increment) {
	if (!debugTop || !enableFlag) {
		return;
	}

	if (increment) {
		iterationCounter++;
	}
	std::ofstream dbg("32bitGBV.log", std::ios::app);
	dbg << functionName << " call #" << iterationCounter;
	outputFunc(dbg);
	dbg << std::endl;
	dbg.close();

	// Check if we should disable all debug flags
	checkAndDisableAllDebugFlags();
}

// Check if all debug iterations have reached the threshold (3) and disable all flags
void checkAndDisableAllDebugFlags() {
	const uint64_t THRESHOLD = 3;

	if (calculateSliceIteration >= THRESHOLD &&
	    getNextSliceIteration >= THRESHOLD &&
	    calculateNodeClipPreciseIteration >= THRESHOLD &&
	    calculateNodeInnerIteration >= THRESHOLD &&
	    getScoreBeforeStartIteration >= THRESHOLD &&
	    mergeTwoSlices2InputIteration >= THRESHOLD &&
	    mergeTwoSlices4InputIteration >= THRESHOLD &&
	    differenceMasksBitTwiddleIteration >= THRESHOLD &&
	    flattenWordSliceIteration >= THRESHOLD) {

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