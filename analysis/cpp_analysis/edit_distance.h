#ifndef RANK_EPISTASIS_EDIT_DISTANCE_H
#define RANK_EPISTASIS_EDIT_DISTANCE_H

template <typename T>
class Matrix2D{
private:
    std::vector<T> vec;
    size_t num_rows;
    size_t num_cols;
public:
    Matrix2D(size_t r, size_t c, T val){
        num_rows = r;
        num_cols = c;
        vec.resize(num_rows * num_cols, val);
    }
    inline T Get(size_t row, size_t col){
        if(row >= num_rows || col >= num_cols){
            std::cerr << "Error! Matrix index out of bounds!" << std::endl;
            std::cerr << "Requested: (" << row << ", " << col << ")" << std::endl;
            std::cerr << "Matrix size: (" << num_rows << ", " << num_cols << ")" << std::endl;
        }
        return vec[row * num_cols + col];
    }
    void Set(int row, int col, T val){
        vec[row * num_cols + col] = val;
    }
    void Print(std::ostream& os = std::cout){
        for(size_t row = 0; row < num_rows; ++row){
            for(size_t col = 0; col < num_cols; ++col){
                std::cout << Get(row, col) << " ";
            }
            std::cout << std::endl;
        }
    }
    void PrintWithLabels(const std::vector<T>& row_labels, const std::vector<T>& col_labels,
            std::ostream& os =  std::cout){
        std::cout << "xxx xxx xxx ";
        for(size_t col = 0; col < num_cols - 2; ++col){
            if(col_labels[col] < 10) std::cout << "  " << col_labels[col] << " ";
            else if(col_labels[col] < 100) std::cout << " " << col_labels[col] << " ";
            else std::cout << col_labels[col] << " ";
        }
        std::cout << std::endl;
        for(size_t row = 0; row < num_rows; ++row){
            if(row <= 1) std::cout << "xxx ";
            else{
                if(row_labels[row - 2] < 10) std::cout << "  " << row_labels[row - 2] << " ";
                else if(row_labels[row - 2] < 100) std::cout << " " << row_labels[row - 2] << " ";
                else std::cout <<  row_labels[row - 2] << " ";
            }
            for(size_t col = 0; col < num_cols; ++col){
                if(Get(row, col) < 10) std::cout << "  " << Get(row, col) << " ";
                else if(Get(row, col) < 100) std::cout << " " << Get(row, col) << " ";
                else std::cout << Get(row, col) << " ";
            }
            std::cout << std::endl;
        }
    }
};


enum EditDistanceMetric{
    kLevenshtein = 0,
    kDamerau_Levenshtein = 1
};

// Levenshtein 
// Implemented from psuedocode at https://en.wikipedia.org/wiki/Levenshtein_distance
double EditDistance_L(const std::vector<size_t>& vec_a, const std::vector<size_t>& vec_b){
    Matrix2D<size_t> matrix(vec_a.size() + 1, vec_b.size() + 1, 0);
    for(size_t i = 0; i <= vec_a.size(); ++i) matrix.Set(i, 0, i);
    for(size_t j = 0; j <= vec_b.size(); ++j) matrix.Set(0, j, j);
    size_t sub_cost = 0;
    for(size_t j = 1; j <= vec_b.size(); ++j){
        for(size_t i = 1; i <= vec_a.size(); ++i){
            sub_cost = (vec_a[i-1] == vec_b[j-1]) ? 0 : 1;
            matrix.Set(i, j, std::min({
                matrix.Get(i, j - 1) + 1,
                matrix.Get(i - 1, j) + 1,
                matrix.Get(i - 1, j - 1) + sub_cost }));
        }
    }
    //matrix.Print();
    return matrix.Get(vec_a.size() - 1, vec_b.size() - 1);
}

// Damerau-Levenshtein
// From psuedocode at https://en.wikipedia.org/wiki/Damerau%E2%80%93Levenshtein_distance 
// Also useful for understanding: https://www.lemoda.net/text-fuzzy/damerau-levenshtein/index.html
double EditDistance_DL(const std::vector<size_t>& vec_a, const std::vector<size_t>& vec_b){
    std::vector<size_t> row_labels = {0, 0};
    for(size_t x : vec_a) row_labels.push_back(x);
    std::vector<size_t> col_labels = {0};
    for(size_t x : vec_b) col_labels.push_back(x);
    
    size_t max_a = *std::max_element(vec_a.begin(), vec_a.end());
    size_t max_b = *std::max_element(vec_b.begin(), vec_b.end());
    std::vector<size_t> da(std::max(max_a, max_b) + 1, 1);
    Matrix2D<size_t> matrix(vec_a.size() + 2, vec_b.size() + 2, 0);
    size_t max_dist = vec_a.size() + vec_b.size();    
    matrix.Set(0, 0, max_dist);
    for(size_t i = 0; i <= vec_a.size() + 1; ++i){
        matrix.Set(i, 0, max_dist);
        matrix.Set(i, 1, i - 1);
    }
    for(size_t j = 1; j <= vec_b.size() + 1; ++j){
        matrix.Set(0, j, max_dist);
        matrix.Set(1, j, j - 1);
    }


    for(size_t i = 2; i <= vec_a.size() + 1; ++i){
        size_t db = 1;
        for(size_t j = 2; j <= vec_b.size() + 1; ++j){
            size_t k = da[vec_b[j - 2]]; // Last row we saw current symbol in vec_b
            size_t l = db;
            size_t cost = 0;
            if(vec_a[i - 2] == vec_b[j - 2]){
                cost = 0;
                db = j;
            }
            else
                cost = 1;
            matrix.Set(i, j, std::min({
                matrix.Get(i - 1, j - 1) + cost,  // Substitution
                matrix.Get(i, j - 1) + 1, // Insertion
                matrix.Get(i - 1, j) + 1, // Deletion
                matrix.Get(k - 1, l - 1) + (i - k - 1) + 1 + (j - l - 1) // Transposition
            }));
        }
        da[vec_a[i - 2]] = i; // The last row we current symbol win vec_a is this one!
    }
    //matrix.PrintWithLabels(vec_a, vec_b);
    return matrix.Get(vec_a.size() + 1, vec_b.size() + 1);
}

// General edit distance function that will direct you to the specified metric
double EditDistance(const std::vector<size_t>& vec_a, const std::vector<size_t>& vec_b, 
        EditDistanceMetric metric){
    switch(metric){
        case kLevenshtein:
            return EditDistance_L(vec_a, vec_b);    
        break;
        case kDamerau_Levenshtein:
            return EditDistance_DL(vec_a, vec_b);    
        break;
        default:
            std::cerr << "Error! Unknown edit distance metric!" << std::endl;
            exit(-1);
            return 0;
        break;
    }
    return 0;
}


#endif
