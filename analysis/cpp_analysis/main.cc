// Standard library
#include <iostream>
#include <vector>
#include <algorithm>
// Empirical
#include "config/ArgManager.h"
#include "config/command_line.h"
// Local
#include "./organism.h"
#include "./file_io.h"
#include "./nk.h"
#include "./edit_distance.h"
#include "./config.h"

void RankOrgs(std::vector<Organism>& orgs){
    for(size_t org_idx = 0; org_idx < orgs.size(); ++org_idx){
        orgs[org_idx].SetID(org_idx);
    }
}

void SortOrgs(std::vector<Organism>& orgs){
    std::stable_sort(orgs.begin(), orgs.end(), [](const Organism& org_a, const Organism& org_b){
        return org_a.GetScore() < org_b.GetScore();
    });
}
// Run a few known examples to check output 
void SanityCheck(){
    // To compare against https://en.wikipedia.org/wiki/Levenshtein_distance
    double dist = EditDistance_L({'s','i','t','t','i','n','g'}, {'k','i','t','t','e','n'});
    std::cout << "Distance: " << dist << std::endl;
    dist = EditDistance_L({'s','u','n','d','a','y'}, {'s','a','t','u','r','d', 'a', 'y'});
    std::cout << "Distance: " << dist << std::endl;
    dist = EditDistance_DL({'a', 'n', ' ', 'a', 'c', 't'}, {'a', ' ', 'c', 'a', 't'});
    std::cout << "Distance: " << dist << std::endl;
    dist = EditDistance_DL({'a', ' ', 'c', 'a', 't'}, {'a', 'n', ' ', 'a', 'c', 't'}); 
    std::cout << "Distance: " << dist << std::endl;
    dist = EditDistance_DL({'a', ' ', 'c', 'a', 't'}, {'a', ' ', 'a', 'b', 'c', 't'}); 
    std::cout << "Distance: " << dist << std::endl;
}


int main(int argc, char* argv[]){
    const std::string config_filename = "analysis_config.cfg";
    // Use command line args from Empirical
    auto args = emp::cl::ArgManager(argc, argv);
    AnalysisConfig config;
    config.Read(config_filename);
    if (args.ProcessConfigOptions(config, std::cout, config_filename, "analysis_config-macros.h")
            == false)
        exit(0);
    if (args.TestUnknown() == false) exit(0); // If there are leftover args, throw an error. 
    // Slurp in variables from config file
    const size_t N =                            (size_t)        config.N();
    const size_t K =                            (size_t)        config.K();
    const size_t gen_start =                    (size_t)        config.GEN_START();
    const size_t gen_end   =                    (size_t)        config.GEN_END();
    const std::string output_filename =         (std::string)   config.OUTPUT_FILENAME();
    const std::string input_filename_prefix =   (std::string)   config.INPUT_FILENAME_PREFIX();
    const std::string input_filename_suffix =   (std::string)   config.INPUT_FILENAME_SUFFIX();
    const size_t edit_distance_metric_tmp =     (size_t)        config.EDIT_DISTANCE_METRIC();
    const EditDistanceMetric edit_distance_metric = (EditDistanceMetric)edit_distance_metric_tmp;
    // Write to screen how the experiment is configured
    std::cout << "==============================" << std::endl;
    std::cout << "|    Current configuration   |" << std::endl;
    std::cout << "==============================" << std::endl;
    config.Write(std::cout);
    std::cout << "==============================\n" << std::endl;
 
     
    // Attempt to open the output file
    std::ofstream fp_out;
    fp_out.open(output_filename, std::ios::out);
    if(!fp_out.is_open()){
        std::cerr << "Unable to open output file: " << output_filename << std::endl;
        exit(-1);
    }
    fp_out << "gen,locus,edit_distance,weighted_edit_distance\n";
    
    // Setup desired NK table
    PositionlessNKTable nk_table(K);
    nk_table.SetValue(BinaryVecToInteger({1,1,1}), 1);
    nk_table.SetValue(BinaryVecToInteger({0,0,0}), 1);
    nk_table.SetValue(BinaryVecToInteger({1,0,1}), 2);
    nk_table.SetValue(BinaryVecToInteger({0,1,0}), 2);
    std::cout << "Using the following NK table:" << std::endl;
    nk_table.Print();
    
    std::vector<Organism> orgs;
    std::stringstream ss;
    for(size_t gen = gen_start; gen < gen_end; ++gen){
        // Load snapshot file
        orgs.clear();
        ss.str("");
        ss << input_filename_prefix << gen << input_filename_suffix;
        size_t num_unique_genomes = LoadOrgsFromSnapshot(ss.str(), orgs); 
        
        std::cout 
            << gen << " "
            << "(" << num_unique_genomes << " / " 
            << orgs.size() << " )" 
            << std::endl;
        
        // Score the organisms
        for(size_t org_idx = 0; org_idx < orgs.size(); ++org_idx){
            orgs[org_idx].Score(K, nk_table);
        }
        // Sort the organisms based on score and then label them
        SortOrgs(orgs);
        RankOrgs(orgs);
 
        // Extract ranks 
        std::vector<size_t> rank_vec(orgs.size(), 0);
        std::vector<size_t> rank_vec_mutated(orgs.size(), 0);
        for(size_t org_idx = 0; org_idx < orgs.size(); ++org_idx){
            rank_vec[org_idx] = orgs[org_idx].GetID();
        }

        // Create a genome_weighting_factor, num. unique genomes / num. total genomes
        double genome_weighting_factor = ((double)num_unique_genomes) / orgs.size(); 

        double edit_distance = 0;
        // Create, score, and sort some mutated organisms
        std::vector<Organism> orgs_mutated(orgs.size());
        for(size_t locus_idx = 0; locus_idx < N; ++locus_idx){ 
            for(size_t org_idx = 0; org_idx < orgs.size(); ++org_idx){
                orgs_mutated[org_idx] = orgs[org_idx].GetMutatedCopy(locus_idx);
                orgs_mutated[org_idx].Score(K, nk_table);
            }
            SortOrgs(orgs_mutated);
            for(size_t org_idx = 0; org_idx < orgs.size(); ++org_idx){
                rank_vec_mutated[org_idx] = orgs_mutated[org_idx].GetID();
            }
            //std::cout << "mutants (" << gen << "," << locus_idx << ") " << std::endl;
            //for (int i = 0; i < orgs.size(); i++) {
            //    std::cout << "[" << rank_vec_mutated[i] << "]"
            //        << orgs_mutated[i].GetScore() / N << " ";
            //}        
            //std::cout << std::endl;
            edit_distance = EditDistance(rank_vec, rank_vec_mutated, edit_distance_metric);
            fp_out << gen << ","
                << locus_idx << ","
                << edit_distance << ","
                << edit_distance * genome_weighting_factor << "\n";
            //std::cout << "edit distance (" << gen << "," << locus_idx << ") " << edit_distance
            //    << std::endl; 
        }
    }
    fp_out.close();
    return 0;
}
