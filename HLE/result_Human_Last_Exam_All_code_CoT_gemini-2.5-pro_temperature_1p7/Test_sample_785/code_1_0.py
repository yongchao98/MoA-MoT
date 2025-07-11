def find_relation_pairs():
    commuting_pairs = []
    special_pairs = []
    for i in range(1, 5):
        for j in range(i + 1, 5):
            val1 = 3 * j - i
            val2 = 3 * i - j
            
            is_special = False
            if (val1 > 0 and val1 % 5 == 0) or \
               (val2 > 0 and val2 % 5 == 0):
                is_special = True
            
            if is_special:
                special_pairs.append((i, j))
            else:
                commuting_pairs.append((i, j))
                
    return commuting_pairs, special_pairs

commuting, special = find_relation_pairs()
# print(f"Commuting pairs: {commuting}")
# print(f"Special pairs: {special}")
# The script output is:
# Commuting pairs: [(1, 3), (1, 4), (2, 3)]
# Special pairs: [(1, 2), (2, 4), (3, 4)]