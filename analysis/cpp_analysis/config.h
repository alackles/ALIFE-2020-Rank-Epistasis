#ifndef RANK_EPISTASIS_ANALYSIS_CONFIG_H
#define RANK_EPISTASIS_ANALYSIS_CONFIG_H

#include "config/config.h"
#include <string>

EMP_BUILD_CONFIG(AnalysisConfig,
    // General Group 
    GROUP(GENERAL, "General settings"),
    VALUE(N,            size_t,     200,    "Number of bits in each genome"),
    VALUE(K,            size_t,     3,      "Number of loci to use for each NK lookup"),
    VALUE(GEN_START,    size_t,     0,      "Minimum generation to analyze (included)"),
    VALUE(GEN_END,      size_t,     500,    "Maximum generation to analye (excluded)"),
    VALUE(OUTPUT_FILENAME,          std::string, "edit_distance.csv", "Path to save output file"),
    VALUE(INPUT_FILENAME_PREFIX,    std::string, "./",  "Input filepath up to generation number"),
    VALUE(INPUT_FILENAME_SUFFIX,    std::string, ".csv","Input filepath after generation number"),
    VALUE(EDIT_DISTANCE_METRIC,     size_t, 0,   "0 for Levenshtein, 1 for Damerau-Levenshtein")
)
#endif
