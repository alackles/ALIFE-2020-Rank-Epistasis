#ifndef RANK_EPISTASIS_NK_H 
#define RANK_EPISTASIS_NK_H 

#include <vector>
#include <iostream>
#include <iomanip>

size_t BinaryVecToInteger(std::vector<unsigned short> vec){
    size_t return_val = 0;
    for(size_t i = 0; i < vec.size(); ++i){
        return_val += (1 <<  i) * vec[i];
    }
    return return_val;
}
void PrintBinary(size_t val, size_t padded_size, std::ostream& os){
    size_t num_digits = 0; 
    while(val >> num_digits != 0){
        ++num_digits;
    }
    for(size_t i = num_digits; i < padded_size; ++i){
        std::cout << "0"; 
    }
    for(size_t i = 0; i < num_digits; ++i){
        std::cout << ((val >> (num_digits - i - 1)) & 1);  
    }
}

class PositionlessNKTable{
private:
    std::vector<double> value_vec;
    size_t K;
public:
    PositionlessNKTable(size_t k, double val = 0){
        value_vec.resize(2 << (k - 1), val);
        K = k;
    }
    void SetValue(size_t idx, double val){
        if(idx < 0 || idx >= value_vec.size()){
            std::cerr << "Error! Attempted to set out-of-bounds index in NKTable!" 
                      << " Index: " << idx  << "!" << std::endl;
            exit(-1);
        }
        value_vec[idx] = val;
    }
    size_t GetValue(size_t idx){
        if(idx < 0 || idx >= value_vec.size()){
            std::cerr << "Error! Attempted to get out-of-bounds index in NKTable!" 
                      << " Index: " << idx  << "!" << std::endl;
            exit(-1);
        }
        return value_vec[idx];
    }
    void Print(std::ostream& os = std::cout){
        for(size_t i = 0; i < value_vec.size(); ++i){
            PrintBinary(i, K, os);
            std::cout << std::setbase(10) << " -> " << value_vec[i] << std::endl;
        } 
    }
};


#endif
