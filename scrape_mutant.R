# Clear!
rm(list = ls())

# Libraries
library(dplyr)
library(hash)

# Config vars
file_of_filenames = 'file_lists/files_TRD_0_final.txt'
filename_list_prefix = strsplit(file_of_filenames, '.txt')[[1]][1]
filename_seed_prefix = '/'
filename_seed_suffix = '/mutant_fitness.csv'
filename_out_prefix = 'data/parts/mutant_fitness_TRD_0_final__'
min_seed = 501
max_seed = 550

# Load in the list of filenames
filename_vec = as.character(read.csv(file_of_filenames, header = F)[,1])
cat(paste0('Found ', length(filename_vec), ' files to scrape!\n'))

# Grab the appropriate column names
column_names = c(
    'update',
    'org_idx', 
    'num_mutations',
    'fitness_avg',
    'fitness_max',
    'fitness_min',
    'condition', 'seed','FIX', 'K', 'ED', 'MUT', 'TRD', 'VEL')
print('Column names to scrape: ')
print(column_names)

condition_num = 1
# Actually start scraping data
for(filename in filename_vec){
    cat(paste0(condition_num, ' / ', length(filename_vec), ': ', filename, '\n'))
    data = data.frame(data = matrix(nrow = 0, ncol = length(column_names)))
    colnames(data) = column_names
    cur_vals = hash()
    cur_vals[['condition']] = strsplit(filename, '__')[[1]][1]
    param_pairs = strsplit(filename, '__')[[1]]
    for(param_pair in param_pairs){
        param_parts = strsplit(param_pair, '_')[[1]]
        if(length(param_parts) > 1){
            cur_vals[[param_parts[1]]] = param_parts[2]
        }
    }
    if(is.null(cur_vals[['TRD']])){
        cur_vals[['TRD']] = 0
        cur_vals[['VEL']] = 0
    }
    for(seed in (min_seed + 1):max_seed){
        filename_seed = paste0(
                filename,
                filename_seed_prefix,
                seed,
                filename_seed_suffix)
        cat(paste0(seed, ' '))
        if(!file.exists(filename_seed)){
            cat('\nMissing file: ', filename_seed, '\n')
            next
        }
        data_seed = cbind(
            read.csv(filename_seed), 
            cur_vals[['condition']],
            seed, 
            cur_vals[['FIX']],
            cur_vals[['K']],
            cur_vals[['ED']],
            cur_vals[['MUT']],
            cur_vals[['TRD']],
            cur_vals[['VEL']]
        )
        data = rbind(data, data_seed)
    }
    cat('\n')
    colnames(data) = column_names
    write.csv(data, file = paste0(filename_out_prefix, condition_num, '.csv'))
    condition_num = condition_num + 1
}
cat(paste0('Done! Saved results to: ', filename_out_prefix, 'n.csv\n'))
