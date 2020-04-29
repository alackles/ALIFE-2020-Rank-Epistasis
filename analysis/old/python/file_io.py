def load_orgs_from_snapshot(filename):
    org_L = []
    with open(filename, 'r') as fp:
        fp.readline()
        for line in fp:
            line = line.strip()
            if line == '':
                continue
            line_part_L = line.split(',')
            genome_size = int(line_part_L[0])
            s = ''
            for offset in range(genome_size):
                s += line_part_L[1 + offset]
            org_L.append(s.strip('"'))
    return org_L

if __name__ == '__main__': 
    org_L = load_orgs_from_snapshot('../MABE/snapshots/flat/snapshot_organisms_1.csv')
    for x in org_L:
        print(x)
