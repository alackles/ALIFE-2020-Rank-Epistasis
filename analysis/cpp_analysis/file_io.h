#ifndef RANK_EPISTASIS_FILE_IO_h
#define RANK_EPISTASIS_FILE_IO_h

#include <fstream>
#include <string>
#include <unordered_set>
#include <sstream>

#include "organism.h"

// Load in the specified snapshot (should be a .csv)
//  and add each organism to the passed vector
// Returns the number of UNIQUE genomes loaded
size_t LoadOrgsFromSnapshot(std::string filename, std::vector<Organism>& org_vec){
    std::ifstream fp;
    std::string line;
    fp.open(filename, std::ios::in);
    // If we failed to open file, error out
    if(!fp.is_open()){
        std::cerr << "Error! Unable to load snapshot file: " << filename << std::endl;
        exit(-1);
    }
    
    bool in_genome = false;
    bool is_organism = false;
    size_t org_idx = 0;
    std::vector<unsigned short> cur_genome;
    std::unordered_set<std::string> genome_set;
    std::stringstream ss;
    // Iterate through each line
    while(getline(fp,line)){
        if(line == "" || line == "\n")
            continue;
        in_genome = false;
        is_organism = false;
        for(size_t i = 0; i < line.length(); ++i){
            if(line[i] == '"'){
                if(!in_genome){
                    org_vec.push_back(Organism());
                    cur_genome.clear(); 
                    ss.str("");
                }
                in_genome = !in_genome;
                is_organism = true; 
            }
            else if(in_genome){
                ss << line[i];
                if(line[i] == '1'){
                    org_vec[org_idx].PushGene(1);
                    cur_genome.push_back(1);
                }
                else if(line[i] == '0'){
                    org_vec[org_idx].PushGene(0);
                    cur_genome.push_back(0);
                }
            }
        } 
        if(is_organism){       
            ++org_idx;
            genome_set.insert(ss.str());
        }
    }
    fp.close();
    return genome_set.size();
}

#endif
