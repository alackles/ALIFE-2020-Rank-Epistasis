snapshot_dir = './archivist_data/'
num_gens = 200
num_traits = 9 

with open('org_data.csv', 'w') as out_fp:
    out_fp.write('gen,id')
    for trait_id in range(num_traits):
        out_fp.write(',trait_' + str(trait_id))
    out_fp.write(',int_val\n')
  
    for cur_gen in range(num_gens + 1):
        with open(snapshot_dir + 'snapshot_organisms_' + str(cur_gen) + '.csv', 'r') as in_fp:
            for line in in_fp:
                line = line.strip()
                parts_L = line.split('"')
                if len(parts_L) != 3:
                    continue
                trait_L = parts_L[1].split(',')
                cur_id = int(parts_L[-1].split(',')[-1])
                out_fp.write(str(cur_gen) + ',' + str(cur_id))
                trait_int_val = 0
                for trait_id in range(num_traits):
                    out_fp.write(',' + trait_L[trait_id])
                    if trait_L[trait_id] == '1':
                        trait_int_val += 1 << (num_traits - 1 - trait_id)
                out_fp.write(',' + str(trait_int_val) + '\n')
