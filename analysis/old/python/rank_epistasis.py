def get_score(bitstring, n, k, nk_table):
    str_copy = bitstring
    score = 0
    for cur_n in range(n): 
        k_str = str_copy[0:k]
        k_score = nk_table[k_str]
        str_copy = str_copy[1:] + str_copy[0]
        score += k_score
    return score

# Get the edit distance between two sequences, a and b
def get_distance(a, b):
    M = []
    for i in range(len(a) + 1):
        M.append([0] * (len(b) + 1))
    for i in range(len(a) + 1):
        M[i][0] = i
    for j in range(len(b) + 1):
        M[0][j] = j
    for j in range(1, len(b) + 1):
        for i in range(1, len(a) + 1):
            sub_cost = 0
            if a[i-1] != b[j-1]: 
                sub_cost = 1
            M[i][j] = min( \
                    M[i - 1][j] + 1, \
                    M[i][j - 1] + 1, \
                    M[i - 1][j - 1] + sub_cost)
    for m in M:
        print(m)
    return M[len(a) - 1][len(b) - 1]
