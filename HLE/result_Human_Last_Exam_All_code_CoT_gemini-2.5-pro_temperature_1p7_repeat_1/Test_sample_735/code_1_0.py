import collections

def is_maximal_product_free(group_name, mult_table, s_indices):
    """
    Verifies if a given set S is a maximal product-free set in a group G.

    Args:
        group_name (str): The name of the group for printing.
        mult_table (list of lists): The multiplication table of the group.
        s_indices (set): A set of integers representing the elements of S.

    Returns:
        bool: True if S is a maximal product-free set, False otherwise.
    """
    num_elements = len(mult_table)
    elements = set(range(num_elements))
    S = set(s_indices)

    # 1. Verify S is product-free: S * S should not intersect S.
    for s1 in S:
        for s2 in S:
            product = mult_table[s1][s2]
            if product in S:
                print(f"FAILED [{group_name}]: Set S={S} is NOT product-free. {s1}*{s2}={product} is in S.")
                return False

    # 2. Verify S is maximal: For any g in G\S, the set S U {g} is NOT product-free.
    g_elements = elements - S
    for g in g_elements:
        s_union_g = S.union({g})
        is_still_product_free = True
        # Check if S U {g} is product-free.
        for x in s_union_g:
            for y in s_union_g:
                product = mult_table[x][y]
                if product in s_union_g:
                    is_still_product_free = False
                    break
            if not is_still_product_free:
                break
        
        # If we found no product in S U {g}, it means S U {g} is still product-free,
        # which violates maximality.
        if is_still_product_free:
            print(f"FAILED [{group_name}]: Set S={S} is NOT maximal. It can be extended with {g}.")
            return False
            
    print(f"SUCCESS [{group_name}]: Set S={S} is a maximal product-free set of size 3.")
    return True

def define_groups_and_verify():
    """
    Defines the 6 known groups and their candidate sets, then runs verification.
    """
    verified_groups = collections.OrderedDict()

    # Group 1: Z_6 (cyclic group of order 6)
    # S = {1, 3, 5}
    z6_table = [[(i + j) % 6 for j in range(6)] for i in range(6)]
    verified_groups['Z_6'] = (z6_table, {1, 3, 5})
    
    # Group 2: S_3 (Dihedral group D_6, order 6)
    # e, r, r^2, s, sr, sr^2 -> 0, 1, 2, 3, 4, 5
    # S = {s, sr, sr^2} -> {3, 4, 5}
    s3_table = [
        [0, 1, 2, 3, 4, 5],
        [1, 2, 0, 5, 3, 4],
        [2, 0, 1, 4, 5, 3],
        [3, 4, 5, 0, 1, 2],
        [4, 5, 3, 2, 0, 1],
        [5, 3, 4, 1, 2, 0]
    ]
    verified_groups['S_3'] = (s3_table, {3, 4, 5})

    # Group 3: Z_8 (cyclic group of order 8)
    # S = {3, 4, 5}
    z8_table = [[(i + j) % 8 for j in range(8)] for i in range(8)]
    verified_groups['Z_8'] = (z8_table, {3, 4, 5})

    # Group 4: D_8 (dihedral group of order 8)
    # e, r, r^2, r^3, s, sr, sr^2, sr^3 -> 0, 1, 2, 3, 4, 5, 6, 7
    # S = {s, sr, r^2} -> {4, 5, 2}
    d8_table = [
        [0, 1, 2, 3, 4, 5, 6, 7],
        [1, 2, 3, 0, 7, 4, 5, 6],
        [2, 3, 0, 1, 6, 7, 4, 5],
        [3, 0, 1, 2, 5, 6, 7, 4],
        [4, 5, 6, 7, 0, 1, 2, 3],
        [5, 6, 7, 4, 3, 0, 1, 2],
        [6, 7, 4, 5, 2, 3, 0, 1],
        [7, 4, 5, 6, 1, 2, 3, 0]
    ]
    verified_groups['D_8'] = (d8_table, {2, 4, 5})

    # Group 5: Z_9 (cyclic group of order 9)
    # S = {2, 3, 7}
    z9_table = [[(i + j) % 9 for j in range(9)] for i in range(9)]
    verified_groups['Z_9'] = (z9_table, {2, 3, 7})

    # Group 6: Z_3 x Z_3
    # Elements (i,j) map to 3*i+j. So 0..8
    # S = {(0,1), (1,0), (2,2)} -> {1, 3, 8}
    z3z3_table = [[0]*9 for _ in range(9)]
    for i in range(9):
        for j in range(9):
            r1, c1 = i // 3, i % 3
            r2, c2 = j // 3, j % 3
            rf, cf = (r1 + r2) % 3, (c1 + c2) % 3
            z3z3_table[i][j] = 3 * rf + cf
    verified_groups['Z_3 x Z_3'] = (z3z3_table, {1, 3, 8})

    # --- Verification ---
    count = 0
    for name, (table, s_set) in verified_groups.items():
        if is_maximal_product_free(name, table, s_set):
            count += 1
            
    print("\n-------------------------------------------------")
    print(f"Final Count: The number of finite groups found is {count}.")
    
if __name__ == '__main__':
    define_groups_and_verify()