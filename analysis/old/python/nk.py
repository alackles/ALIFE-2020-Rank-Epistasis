def get_nk_table(n, k, val = 0):
    nk_table = {}
    for i in range(2 ** k):
        tmp_str = str(bin(i))[2:]
        tmp_str = ('0' * (k - len(tmp_str))) + tmp_str
        nk_table[tmp_str] = val
    return nk_table
