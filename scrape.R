# Clear!
rm(list = ls())

# Libraries
library(dplyr)
library(hash)

# Config vars
file_of_filenames = 'dirs_to_scrape.txt'
filename_list_prefix = strsplit(file_of_filenames, '.txt')[[1]][1]
filename_seed_prefix = '/'
filename_seed_suffix = '/edit_distance.csv'
filename_out = 'edit_distance_scraped.csv'
min_seed = 101
max_seed = 150

# Load in the list of filenames
filename_vec = as.character(read.csv(file_of_filenames, header = F)[,1])
cat(paste0('Found ', length(filename_vec), ' files to scrape!\n'))

# Grab the appropriate column names
column_names = c('condition','update','bit_idx')
# Grab condition vars from filename
example_filename = filename_vec[1]
param_pairs = strsplit(example_filename, '__')[[1]]
for(param_pair in param_pairs){
    param_parts = strsplit(param_pair, '_')[[1]]
    if(length(param_parts) > 1){
        column_names = c(column_names, param_parts[1])
    }
}
column_names = c(column_names, 'edit_distance_mean', 'edit_distance_sd')
column_names = c("update","bit_idx","edit_distance","seed","FIX","K","ED")

# Release vars from memory
rm(example_filename)

# Prepare the data frame to hold the scraped data
data = data.frame(data = matrix(nrow = 0, ncol = length(column_names)))
colnames(data) = column_names
print('Column names to scrape: ')
print(colnames(data))

condition_num = 1
# Actually start scraping data
for(filename in filename_vec){
    cat(paste0(condition_num, ' / ', length(filename_vec), ': ', filename, '\n'))
    condition_num = condition_num + 1
    data_condition = read.csv(paste0(
        filename,
        filename_seed_prefix,
        min_seed,
        filename_seed_suffix
    ))
    cur_vals = hash()
    cur_vals[['condition']] = strsplit(filename, '__')[[1]][1]
    param_pairs = strsplit(filename, '__')[[1]]
    for(param_pair in param_pairs){
        param_parts = strsplit(param_pair, '_')[[1]]
        if(length(param_parts) > 1){
            cur_vals[[param_parts[1]]] = param_parts[2]
        }
    }
    data_condition$seed = min_seed
    if(min_seed != max_seed){
        for(seed in (min_seed + 1):max_seed){
            cat(paste0(seed, ' '))
            data_seed = cbind(
                read.csv(paste0(
                    filename,
                    filename_seed_prefix,
                    min_seed,
                    filename_seed_suffix
                )), seed, 
                    cur_vals[['FIX']],
                    cur_vals[['K']],
                    cur_vals[['ED']]
                    )
            #data_condition = rbind(data_condition, data_seed)
            #data_summary$condition = cur_vals[['condition']]
            #data_summary$FIX = cur_vals[['FIX']]
            #data_summary$K= cur_vals[['K']]
            #data_summary$ED= cur_vals[['ED']]
            data = rbind(data, data_seed)
            #data_condition[is.na(data_condition$seed),'seed'] = seed
        }
    }
    #cat('\n') 
    #data_grouped = dplyr::group_by(data_condition, bit_idx, update)
    #data_summary = dplyr::summarize(data_grouped, 
    #    edit_distance_mean = mean(edit_distance), 
    #    edit_distance_sd = sd(edit_distance))
    #print(colnames(data_summary))
    #cat('\n')
    #data_summary$condition = cur_vals[['condition']]
    #data_summary$FIX = cur_vals[['FIX']]
    #data_summary$K= cur_vals[['K']]
    #data_summary$ED= cur_vals[['ED']]
    #data = rbind(data_summary, data)
    #for(column_name in colnames(data)){
    #    data_summary = cbind(data_summary, cur_vals[[column_name]])
    #}
    #print('data')
    #print(colnames(data))
    #print('data_summary')
    #print(colnames(data_summary))

    #data = rbind(data, data_summary)
    #print(head(data))
    #colnames(data) = column_names
}
write.csv(data, filename_out)
cat(paste0('Done! Saved result to: ', filename_out, '\n'))
