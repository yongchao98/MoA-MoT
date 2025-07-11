import itertools

def is_product_free(s, cayley_table):
    """
    Checks if a subset s (represented as a set of indices) is product-free
    for a group defined by its Cayley table.
    """
    if not s:
        return True
    for x in s:
        for y in s:
            product = cayley_table[x][y]
            if product in s:
                return False
    return True

def find_groups_with_maximal_product_free_sets_of_size_2():
    """
    Finds and counts finite groups that have a maximal by inclusion
    product-free set of size 2.
    """
    groups = {}

    # Define Cayley tables for small groups
    # Identity element is always 0
    # Cyclic group C_n
    groups['C2'] = [[(i + j) % 2 for j in range(2)] for i in range(2)]
    groups['C3'] = [[(i + j) % 3 for j in range(3)] for i in range(3)]
    groups['C4'] = [[(i + j) % 4 for j in range(4)] for i in range(4)]
    groups['C5'] = [[(i + j) % 5 for j in range(5)] for i in range(5)]
    groups['C6'] = [[(i + j) % 6 for j in range(6)] for i in range(6)]
    groups['C7'] = [[(i + j) % 7 for j in range(7)] for i in range(7)]
    
    # Klein four-group (C2 x C2)
    groups['V4'] = [[0, 1, 2, 3], [1, 0, 3, 2], [2, 3, 0, 1], [3, 2, 1, 0]]
    
    # Symmetric group S3 (order 6)
    # e=0, (12)=1, (13)=2, (23)=3, (123)=4, (132)=5
    groups['S3'] = [
        [0, 1, 2, 3, 4, 5], [1, 0, 5, 4, 3, 2], [2, 4, 0, 5, 1, 3],
        [3, 5, 4, 0, 2, 1], [4, 2, 3, 1, 5, 0], [5, 3, 1, 2, 0, 4]
    ]
    
    # Quaternion group Q8 (order 8)
    # 1=0, -1=1, i=2, -i=3, j=4, -j=5, k=6, -k=7
    groups['Q8'] = [
        [0, 1, 2, 3, 4, 5, 6, 7], [1, 0, 3, 2, 5, 4, 7, 6],
        [2, 3, 1, 0, 6, 7, 5, 4], [3, 2, 0, 1, 7, 6, 4, 5],
        [4, 5, 7, 6, 1, 0, 2, 3], [5, 4, 6, 7, 0, 1, 3, 2],
        [6, 7, 4, 5, 3, 2, 1, 0], [7, 6, 5, 4, 2, 3, 0, 1]
    ]

    # Dihedral group D10 (order 10, symmetries of a pentagon)
    # e=0, r=1..4, s=5, sr=6..9
    D10_table = [[0]*10 for _ in range(10)]
    for i in range(10):
        for j in range(10):
            i_is_s, i_pow = (i >= 5), i % 5
            j_is_s, j_pow = (j >= 5), j % 5
            if not i_is_s and not j_is_s: res_pow, res_is_s = (i_pow + j_pow) % 5, False
            elif not i_is_s and j_is_s:  res_pow, res_is_s = (j_pow - i_pow) % 5, True
            elif i_is_s and not j_is_s:  res_pow, res_is_s = (i_pow + j_pow) % 5, True
            else: # i_is_s and j_is_s
                res_pow, res_is_s = (j_pow - i_pow) % 5, False
            D10_table[i][j] = res_pow + 5 * res_is_s
    groups['D10'] = D10_table

    found_groups = []
    
    for name, cayley_table in groups.items():
        n = len(cayley_table)
        identity = 0
        non_identity_elements = [i for i in range(n) if i != identity]
        
        if len(non_identity_elements) < 2:
            continue

        group_has_such_set = False
        # Iterate over all 2-element subsets of non-identity elements
        for s_tuple in itertools.combinations(non_identity_elements, 2):
            s = set(s_tuple)
            
            # 1. Check if S is product-free
            if not is_product_free(s, cayley_table):
                continue
            
            # 2. If it is, check if it's maximal
            is_maximal = True
            elements_to_check = [g for g in range(n) if g not in s]
            for g in elements_to_check:
                s_union_g = s.union({g})
                # If S U {g} is still product-free, then S is not maximal
                if is_product_free(s_union_g, cayley_table):
                    is_maximal = False
                    break
            
            if is_maximal:
                group_has_such_set = True
                break # Found one for this group, move to the next group
        
        if group_has_such_set:
            found_groups.append(name)
    
    print(f"Found {len(found_groups)} groups with the specified property:")
    print(sorted(found_groups))
    print("\nTotal number of groups found:")
    print(len(found_groups))

if __name__ == '__main__':
    find_groups_with_maximal_product_free_sets_of_size_2()