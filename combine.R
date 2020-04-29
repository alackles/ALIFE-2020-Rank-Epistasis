rm(list = ls())

filename_in_prefix = './data/parts/edit_distance_TRD_1_final__'
filename_in_suffix = '.csv'
filename_out = './data/edit_distance_TRD_1_final.csv'
idx_start = 1
idx_stop = 27

data_initialized = F
for(idx in idx_start:idx_stop){
    cat(paste0(idx, ' '))
    filename = paste0(filename_in_prefix, idx, filename_in_suffix)
    if(!file.exists(filename)){
        cat('\nMissing file: ', filename, '\n')
        next
    }
    data_idx = read.csv(filename)
    if(!data_initialized){
        data = data_idx
        data_initialized = T
    } else {
        data = rbind(data, data_idx)
    }
}
cat('\n')
write.csv(data, filename_out)
print(paste0('Done! File saved to: ', filename_out, '!'))
