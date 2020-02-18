N = 200
K = 6

val_map = {}
val_map[int('0b111111', 2)] = '1'
val_map[int('0b000000', 2)] = '1'
val_map[int('0b100001', 2)] = '2'
val_map[int('0b110011', 2)] = '1.66'
val_map[int('0b111011', 2)] = '1.33'
val_map[int('0b110111', 2)] = '1.33'
val_map[int('0b011110', 2)] = '2'
val_map[int('0b001100', 2)] = '1.66'
val_map[int('0b001000', 2)] = '1.33'
val_map[int('0b000100', 2)] = '1.33'

with open('output_table.dat', 'w') as out_fp:
    for n in range(N):
        for k in range(2**K):
            if k in val_map.keys():
                out_fp.write(val_map[k] + ' ')
            else:
                out_fp.write('0 ')
        out_fp.write('\n') 
