#ifndef RANK_EPISTASIS_ORGANISM_H
#define RANK_EPISTASIS_ORGANISM_H

// Standard library
#include <vector>
#include <iostream>
// Local
#include "./nk.h"

class Organism{
private:
    std::vector<unsigned short> genome;
    double score;
    int id;
    bool score_needs_calculated = true;
public:
    Organism():
    score(0), id(-1)
    {
    }
    Organism(const Organism& other){
        score = other.score;
        id = other.id;
        score_needs_calculated = other.score_needs_calculated;
        genome = other.genome;
    }
    void PushGene(unsigned short gene_val){
        genome.push_back(gene_val);
        score_needs_calculated = true;
    }
    void Print(std::ostream& os = std::cout) const{
        os << "[" << id <<  "]<";
        for(size_t i = 0; i < genome.size(); ++i){
            os << genome[i];
        }
        os << ">" << std::endl;
    }
    void Score(size_t K, PositionlessNKTable& nk_table){
        score = 0;
        std::vector<unsigned short> gene;
        gene.resize(K, 0);
        for(size_t idx = 0; idx < genome.size(); ++idx){
            for(size_t k = 0; k < K; ++k){
                gene[k] = genome[(idx + k) % genome.size()];
            }
            score += nk_table.GetValue(BinaryVecToInteger(gene));  
        }   
        score_needs_calculated = false;
    }
    double GetScore() const{
        if(score_needs_calculated){
            std::cerr << "Error! Tried to fetch score before it was calculated!" << std::endl;
            exit(-1);
        }
        return score;
    }
    void SetID(int org_id){
        id = org_id;
    }
    int GetID() const{
        return(id);
    }
    Organism GetMutatedCopy(size_t gene_idx) const{
        Organism org = Organism(*this);
        org.genome[gene_idx] = (org.genome[gene_idx] + 1) & 1;
        org.score_needs_calculated = true;
        return org;
    }
    size_t GetGenomeLength() const{
        return genome.size();
    }
};
#endif
