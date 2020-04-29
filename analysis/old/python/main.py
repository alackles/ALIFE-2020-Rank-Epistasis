from rank_epistasis import *
from file_io import *
from viz import *
from nk import *

#print(get_distance('sitting', 'kitten'))
#print(get_distance('sunday', 'saturday'))
#exit(0)

num_gens = 0
with open('edit_distance.csv', 'w') as out_fp:
    out_fp.write('gen,locus,edit_distance\n')
    for gen in range(0, num_gens + 1):
        print(gen, '/', num_gens)
        org_L = load_orgs_from_snapshot(\
            '../MABE/snapshots/1/flat/snapshot_organisms_' + str(gen) + '.csv')

        if len(org_L) <= 0:
            print('Empty organism list!')
            exit(-1)

        n = len(org_L[0])
        k = 3

        nk_table = get_nk_table(n,k,0)
        nk_table['111'] = 1
        nk_table['000'] = 1
        nk_table['101'] = 2
        nk_table['010'] = 2

        score_L = [0] * len(org_L)
        for i in range(len(org_L)):
            score = get_score(org_L[i], n, k, nk_table)
            score_L[i] = score

        #sorted_orgs = [x for _,x in sorted(zip(score_L,org_L))]
        id_L = [x for x in range(len(org_L))]
        sorted_ids = [x for _,x in sorted(zip(score_L, id_L))]

        mut_score_L = [0] * len(org_L)
        for mut_idx in range(1):
            for org_idx in range(len(org_L)):
                mut_org = org_L[org_idx]
                org = org_L[org_idx]
                if mut_org[mut_idx] == '0':
                    mut_org = org[:mut_idx] + '1' + org[mut_idx + 1:]
                else:
                    mut_org = org[:mut_idx] + '0' + org[mut_idx + 1:]
                score = get_score(mut_org, n, k, nk_table)
                mut_score_L[org_idx] = score
            mut_sorted_ids = [x for _,x in sorted(zip(mut_score_L, id_L))]
            edit_dist = get_distance(sorted_ids, mut_sorted_ids)
            out_fp.write(str(gen) + ',' + str(mut_idx) + ',' + str(edit_dist) + '\n')

        jumps = []
        for i in range(len(org_L)):
            old_rank = i
            org_id = sorted_ids[i]
            new_rank = 0
            for j in range(len(org_L)):
                if mut_sorted_ids[j] == org_id:
                    new_rank = j
                    break
            jumps.append(new_rank - old_rank)
        #print(jumps)


        #visualize(sorted_ids, mut_sorted_ids, 1)
