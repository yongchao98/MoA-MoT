from itertools import combinations

def is_product_free(S, group_mult, elements):
    """Checks if a set S is product-free."""
    for s1_name in S:
        for s2_name in S:
            s1 = elements[s1_name]
            s2 = elements[s2_name]
            product = group_mult(s1, s2)
            if product in [elements[s] for s in S]:
                return False
    return True

def is_maximal(S, G_elements, G_identity, group_mult, elements):
    """Checks if a product-free set S is maximal."""
    G_minus_S_names = [name for name in G_elements if name not in S]
    
    for g_name in G_minus_S_names:
        if g_name == G_identity:
            # Adding the identity always breaks product-freeness,
            # as s * e = s (if e is added) or e * e = e
            continue
            
        S_plus_g_names = list(S) + [g_name]
        
        # Check if S u {g} is product-free. If it is, S is not maximal.
        s_plus_g_is_pf = True
        is_broken = False
        for x_name in S_plus_g_names:
            for y_name in S_plus_g_names:
                x = elements[x_name]
                y = elements[y_name]
                product = group_mult(x, y)
                # Check if product is in the new set
                if product in [elements[name] for name in S_plus_g_names]:
                    s_plus_g_is_pf = False
                    is_broken = True
                    break
            if is_broken:
                break
        
        if s_plus_g_is_pf:
            return False # S is not maximal
            
    return True # S is maximal

def solve():
    """
    Checks small groups for maximal product-free sets of size 3.
    """
    found_groups_count = 0
    
    # --- Group 1: S_3 (Symmetric Group on 3 elements) ---
    def s3_mult(p1, p2):
        # Permutations on {0, 1, 2}
        return (p1[p2[0]], p1[p2[1]], p1[p2[2]])

    s3_e   = (0, 1, 2)
    s3_12  = (1, 0, 2)
    s3_13  = (2, 1, 0)
    s3_23  = (0, 2, 1)
    s3_123 = (1, 2, 0)
    s3_132 = (2, 0, 1)
    
    S3_elements = {
        "e": s3_e, "(12)": s3_12, "(13)": s3_13, "(23)": s3_23,
        "(123)": s3_123, "(132)": s3_132
    }
    s3_non_identity = [name for name in S3_elements if name != "e"]

    s3_has_mpfs3 = False
    for s_tuple in combinations(s3_non_identity, 3):
        S = list(s_tuple)
        if is_product_free(S, s3_mult, S3_elements):
            if is_maximal(S, S3_elements.keys(), "e", s3_mult, S3_elements):
                s3_has_mpfs3 = True
                break
    if s3_has_mpfs3:
        found_groups_count += 1

    # --- Group 2: C_3 x C_3 (Abelian group of order 9) ---
    def c3c3_mult(p1, p2):
        return ((p1[0] + p2[0]) % 3, (p1[1] + p2[1]) % 3)

    C3C3_elements = {f"({i},{j})": (i, j) for i in range(3) for j in range(3)}
    c3c3_non_identity = [name for name in C3C3_elements if name != "(0,0)"]

    c3c3_has_mpfs3 = False
    if len(c3c3_non_identity) >= 3:
        for s_tuple in combinations(c3c3_non_identity, 3):
            S = list(s_tuple)
            if is_product_free(S, c3c3_mult, C3C3_elements):
                if is_maximal(S, C3C3_elements.keys(), "(0,0)", c3c3_mult, C3C3_elements):
                    c3c3_has_mpfs3 = True
                    break
    if c3c3_has_mpfs3:
        found_groups_count += 1
    
    print(f"Based on a search of the groups S_3 and C_3 x C_3, the number of groups found to have a maximal product-free set of size 3 is:")
    print(found_groups_count)

solve()
<<<2>>>