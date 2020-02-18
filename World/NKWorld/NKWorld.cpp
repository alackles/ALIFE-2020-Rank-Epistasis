//  MABE is a product of The Hintze Lab @ MSU
//     for general research information:
//         hintzelab.msu.edu
//     for MABE documentation:
//         github.com/Hintzelab/MABE/wiki
//
//  Copyright (c) 2015 Michigan State University. All rights reserved.
//     to view the full license, visit:
//         github.com/Hintzelab/MABE/wiki/License

// Evaluates agents in an NK landscape
// Landscape can be changing periodically (treadmilling) or static

// Note that brain outputs must be 0 or 1; there is a conversion from negative doubles to positive doubles in the landscape
// If brain outputs are only positive, your organisms will have a bad time

#include "NKWorld.h"
#include "../../Utilities/Random.h"
#include <vector>
#include <map>
#include <iostream>
#include <fstream>

#define PI 3.14159265

std::shared_ptr<ParameterLink<int>> NKWorld::nPL =
Parameters::register_parameter("WORLD_NK-n", 4,
        "number of outputs (e.g. traits, loci)");
std::shared_ptr<ParameterLink<int>> NKWorld::kPL =
Parameters::register_parameter("WORLD_NK-k", 2,
        "Number of sites each site interacts with");
std::shared_ptr<ParameterLink<bool>> NKWorld::treadmillPL =
Parameters::register_parameter("WORLD_NK-treadmill", false,
        "whether landscape should treadmill over time. "
        "0 = static landscape, 1 = treadmilling landscape");
std::shared_ptr<ParameterLink<bool>> NKWorld::readNKTablePL =
Parameters::register_parameter("WORLD_NK-readNKTable", false,
        "If true, reads the NK table from the file "
        "specified by inputNKTableFilename "
        "0 = random, 1 = load from file");
std::shared_ptr<ParameterLink<std::string>> NKWorld::inputNKTableFilenamePL =
Parameters::register_parameter("WORLD_NK-inputNKTableFilename", (std::string) "./table_csv",
        "If readNKTable is 1, which file should we use "
        "to load the table?");
std::shared_ptr<ParameterLink<bool>> NKWorld::writeNKTablePL =
Parameters::register_parameter("WORLD_NK-writeNKTable", true,
        "do you want the NK table to be output for each replicate? "
        "0 = no output, 1 = please output");
std::shared_ptr<ParameterLink<double>> NKWorld::velocityPL =
Parameters::register_parameter("WORLD_NK-velocity", 0.01,
        "If treadmilling, how fast should it treadmill? "
        "Smaller values = slower treadmill");
std::shared_ptr<ParameterLink<int>> NKWorld::evaluationsPerGenerationPL =
Parameters::register_parameter("WORLD_NK-evaluationsPerGeneration", 1,
        "Number of times to test each Genome per "
        "generation (useful with non-deterministic "
        "brains)");
std::shared_ptr<ParameterLink<std::string>> NKWorld::groupNamePL =
Parameters::register_parameter("WORLD_NK_NAMES-groupNameSpace",
        (std::string) "root::",
        "namespace of group to be evaluated");
std::shared_ptr<ParameterLink<std::string>> NKWorld::brainNamePL =
Parameters::register_parameter(
        "WORLD_NK_NAMES-brainNameSpace", (std::string) "root::",
        "namespace for parameters used to define brain");

NKWorld::NKWorld(std::shared_ptr<ParametersTable> PT_)
    : AbstractWorld(PT_) {

    // localize N & K parameters
    N = nPL->get(PT);
    K = kPL->get(PT);

    // generate NK lookup table
    // dimensions: N x 2^K
    // each value is a randomly generated pair of doubles, each in [0.0,1.0]
    // represents weighting on the fitness fcn for each k-tuple in the brain output state
    NKTable.clear();
    NKTable.resize(N);
    if(readNKTablePL->get(PT)){
        std::cout << "Attempting to read NK table from file: " << inputNKTableFilenamePL->get(PT) 
                  << std::endl;
        std::ifstream tableFP;
        tableFP.open(inputNKTableFilenamePL->get(PT), std::ios::in);
        if(!tableFP.is_open()){
            std::cerr << "ERROR! Unable to open file: " 
                      << inputNKTableFilenamePL->get(PT) << std::endl;
            exit(-1);
        }
        double cur_val = 0;
        std::cout << "NK table:" << std::endl;
        for(size_t n = 0; n < N; ++n){
            NKTable[n].resize(1<<K);
            for(int k=0;k<(1<<K);k++){
                tableFP >> cur_val;
                NKTable[n][k]= std::pair<double,double>(cur_val, 0);
                if(k != 0) std::cout << " ";
                std::cout << cur_val;
            }
            std::cout << std::endl;
        } 
    }
    else{
        for(int n=0;n<N;n++){
            NKTable[n].resize(1<<K);
            for(int k=0;k<(1<<K);k++){
                NKTable[n][k]= std::pair<double,double>(Random::getDouble(0.0,1.0),Random::getDouble(0.0,1.0));
            }
        }  
    }

    if (writeNKTablePL->get(PT)) {
        std::ofstream NKTable_csv;
        NKTable_csv.open("NKTable.csv");
        for(int k=0;k<(1<<K);k++){
            for(int n=0;n<N;n++){
                NKTable_csv << NKTable[n][k].first;
                // we don't want commas on the last one
                if (n < N-1) {
                    NKTable_csv << ",";
                } else {
                    NKTable_csv << "\n";
                }
            }
        }
        NKTable_csv.close();
    }

    // columns to be added to ave file
    popFileColumns.clear();
    popFileColumns.push_back("score");
    popFileColumns.push_back("score_VAR"); // specifies to also record the
    // variance (performed automatically
    // because _VAR)
}

// create angular sin function for more even fitness distribution
double NKWorld::triangleSin(double x) {
    double Y = 0.0;
    for (int i = 0; i<N; i++) {
        double n = (2*i) + 1;
        Y += pow(-1, (double)i)*pow(n, -2.0)*sin(n*x);
    }
    return (0.25*PI)*Y;
}

void NKWorld::evaluateSolo(std::shared_ptr<Organism> org, int analyze,
        int visualize, int debug) {
    auto brain = org->brains[brainNamePL->get(PT)];
    double t = Global::update*(velocityPL->get(PT));
    for (int r = 0; r < evaluationsPerGenerationPL->get(PT); r++) {

        brain->resetBrain();
        brain->update();

        // fitness function
        double W = 0.0;
        for (int n=0;n<N;n++) {
            //std::cout << "(" << n << ") ";
            int val = 0;
            for (int k=0; k<K; k++) {
                //std::cout <<  (brain->readOutput((n+k)%N) > 0.0);
                val = (val<<1) + (brain->readOutput((n+k)%N) > 0.0); // convert k adjacent sites to integer for indexing NK table
            }
            if (treadmillPL->get(PT)) {
                // formula for localValue generated by Arend Hintze
                double alpha = NKTable[n][val].first;   
                double beta = NKTable[n][val].second;
                double localValue = (1.0 + triangleSin((t*(beta+0.5))+(alpha*PI*2.0)))/2.0;
                W += localValue;
            } else {
                W += NKTable[n][val].first;
            }
            //std::cout << " " << NKTable[n][val].first << " ";
        }
        //std::cout << std::endl;
        double score = W/(double)N;

        org->dataMap.append("score", score);
        if (visualize) {
            std::string filename = "postscore_" + Global::initPopPL->get(PT) + ".csv";
            std::ofstream fileout;
            fileout.open(filename,std::ios::app);
            fileout << org->ID << "," << score << "\n";
            fileout.close();
        }
    }
}

