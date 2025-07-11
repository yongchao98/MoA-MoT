import itertools

def verify_maximal_product_free_set(group_name, element_names, cayley_table, s_indices):
    """
    Verifies if a given set S is a maximal product-free set in group G.

    Args:
        group_name (str): The name of the group for printing.
        element_names (list): A list of strings representing the elements.
        cayley_table (list of lists): The multiplication table of the group.
        s_indices (set): A set of indices corresponding to the product-free set S.

    Returns:
        bool: True if S is a maximal product-free set of size 2, False otherwise.
    """
    print(f"--- Verifying {group_name} ---")
    s_elements = {element_names[i] for i in s_indices}
    print(f"Candidate set S = {s_elements}")

    # 1. Verify S is product-free
    # For all a,b in S, a*b must not be in S.
    for i in s_indices:
        for j in s_indices:
            product_idx = cayley_table[i][j]
            if product_idx in s_indices:
                print(f"FAILED: Set is not product-free.")
                print(f"Reason: {element_names[i]} * {element_names[j]} = {element_names[product_idx]}, which is in S.")
                return False

    print("OK: Set is product-free.")

    # 2. Verify S is maximal
    # For any g not in S, the set S U {g} must NOT be product-free.
    all_indices = set(range(len(element_names)))
    g_indices = all_indices - s_indices
    
    for g_idx in g_indices:
        s_union_g_indices = s_indices.union({g_idx})
        
        # Check if s_union_g_indices is product-free
        is_still_product_free = True
        for i in s_union_g_indices:
            for j in s_union_g_indices:
                product_idx = cayley_table[i][j]
                if product_idx in s_union_g_indices:
                    # Found a product within the new set, so it's NOT product-free.
                    # This is what we want for maximality.
                    is_still_product_free = False
                    break
            if not is_still_product_free:
                break
        
        # If the loop finished without finding a product inside the set,
        # it means S U {g} is still product-free, so S was not maximal.
        if is_still_product_free:
            print(f"FAILED: Set is not maximal.")
            g_element = element_names[g_idx]
            print(f"Reason: S U {{'{g_element}'}} is still a product-free set.")
            return False

    print(f"OK: Set is maximal.")
    print(f"SUCCESS: {group_name} contains a maximal product-free set of size 2.")
    return True

def main():
    verified_groups = []

    # Group 1: C4 (Cyclic group of order 4)
    c4_elements = ['e', 'g', 'g^2', 'g^3']
    c4_table = [[(i + j) % 4 for j in range(4)] for i in range(4)]
    c4_s = {1, 3}  # S = {g, g^3}
    if verify_maximal_product_free_set('C4', c4_elements, c4_table, c4_s):
        verified_groups.append('C4')
        
    # Group 2: C5 (Cyclic group of order 5)
    c5_elements = ['e', 'g', 'g^2', 'g^3', 'g^4']
    c5_table = [[(i + j) % 5 for j in range(5)] for i in range(5)]
    c5_s = {2, 3}  # S = {g^2, g^3}
    if verify_maximal_product_free_set('C5', c5_elements, c5_table, c5_s):
        verified_groups.append('C5')

    # Group 3: C2 x C2 (Klein four-group)
    v4_elements = ['e', 'a', 'b', 'c']
    # e=(0,0), a=(1,0), b=(0,1), c=(1,1)
    v4_map = [(0,0), (1,0), (0,1), (1,1)]
    v4_rev_map = {(0,0):0, (1,0):1, (0,1):2, (1,1):3}
    v4_table = [[0]*4 for _ in range(4)]
    for i in range(4):
        for j in range(4):
            res_tuple = ((v4_map[i][0] + v4_map[j][0])%2, (v4_map[i][1] + v4_map[j][1])%2)
            v4_table[i][j] = v4_rev_map[res_tuple]
    v4_s = {1, 2}  # S = {a, b}
    if verify_maximal_product_free_set('C2 x C2', v4_elements, v4_table, v4_s):
        verified_groups.append('C2 x C2')
        
    # Group 4: S3 (Symmetric group on 3 elements)
    s3_elements = ['e', 'x', 'x^2', 'y', 'xy', 'x^2y']
    # e=0, x=1, x^2=2, y=3, xy=4, x^2y=5
    # Relations: x^3=e, y^2=e, yx = x^2y
    s3_table = [
    #    e,   x, x^2,   y,  xy, x^2y
        [0,   1,   2,   3,   4,   5], # e
        [1,   2,   0,   4,   5,   3], # x
        [2,   0,   1,   5,   3,   4], # x^2
        [3,   5,   4,   0,   2,   1], # y
        [4,   3,   5,   1,   0,   2], # xy
        [5,   4,   3,   2,   1,   0]  # x^2y
    ]
    s3_s = {1, 3}  # S = {x, y}
    if verify_maximal_product_free_set('S3', s3_elements, s3_table, s3_s):
        verified_groups.append('S3')
    
    print("\n--- Summary ---")
    print("The complete list of finite groups containing a maximal product-free set of size 2 is:")
    print(f"{verified_groups}")
    
    equation_parts = [f"1 ({name})" for name in verified_groups]
    total = len(verified_groups)
    print("The final count is:")
    print(f"{' + '.join(equation_parts)} = {total}")

if __name__ == "__main__":
    main()