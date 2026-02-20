#pragma once
#include <cstdint>
#include <functional>
#include <string>

extern bool enableCalculateSliceDebug;
extern uint64_t calculateSliceIteration;

extern bool enableGetNextSliceDebug;
extern uint64_t getNextSliceIteration;

extern bool enableCalculateNodeClipPreciseDebug;
extern uint64_t calculateNodeClipPreciseIteration;

extern bool enableCalculateNodeInnerDebug;
extern uint64_t calculateNodeInnerIteration;

extern bool enableGetScoreBeforeStartDebug;
extern uint64_t getScoreBeforeStartIteration;

extern bool enableMergeTwoSlices2InputDebug;
extern uint64_t mergeTwoSlices2InputIteration;

extern bool enableMergeTwoSlices4InputDebug;
extern uint64_t mergeTwoSlices4InputIteration;

extern bool enableDifferenceMasksBitTwiddleDebug;
extern uint64_t differenceMasksBitTwiddleIteration;

extern bool enableFlattenWordSliceDebug;
extern uint64_t flattenWordSliceIteration;

extern uint64_t EqVectorIteration;


extern bool debugTop;

// Helper function for debug logging
// Usage: DEBUG_LOG("functionName", enableFlag, iterationCounter,
//                  [&](std::ostream& dbg) { dbg << "custom output"; });
// Set increment=false for "Post Op" or follow-up logs that should use the same iteration number
void DEBUG_LOG(const std::string& functionName,
               bool& enableFlag,
               uint64_t& iterationCounter,
               std::function<void(std::ostream&)> outputFunc,
               bool increment = true);

// Check if all debug iterations have reached the threshold and disable all flags
void checkAndDisableAllDebugFlags();
