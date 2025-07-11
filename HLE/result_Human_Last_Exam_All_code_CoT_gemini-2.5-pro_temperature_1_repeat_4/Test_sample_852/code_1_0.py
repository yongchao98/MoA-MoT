import itertools

def get_group_elements(dims):
    """Generates elements of a direct product of cyclic groups Z_n1 x Z_n2 x ..."""
    pools = [range(d) for d in dims]
    return list(itertools.product(*pools))

def add_elements(p1, p2, dims):
    """Adds two group elements component-wise with respective moduli."""
    return tuple((p1[i] + p2[i]) % dims[i] for i in range(len(dims)))

def is_sum_free(S, G, dims):
    """Checks if a set S is sum-free."""
    s_set = set(S)
    for s1 in S:
        for s2 in S:
            s_sum = add_elements(s1, s2, dims)
            if s_sum in s_set:
                return False
    return True

def is_maximal(S, G, dims):
    """Checks if a sum-free set S is maximal."""
    s_set = set(S)
    g_minus_s = [g for g in G if g not in s_set]
    
    for g in g_minus_s:
        # Check if S U {g} is sum-free. If it is for any g, then S is not maximal.
        new_set_list = S + [g]
        if is_sum_free(new_set_list, G, dims):
            return False
    return True

def find_k_S(S, G, dims):
    """Finds the set k(S) = {g in G | 2g in S}."""
    s_set = set(S)
    k_S = []
    for g in G:
        g_doubled = add_elements(g, g, dims)
        if g_doubled in s_set:
            k_S.append(g)
    return k_S

def main():
    """
    Main function to find the smallest abelian group and set S satisfying the condition.
    This function demonstrates the solution for G = Z_2 x Z_2 x Z_4.
    """
    group_dims = (2, 2, 4)
    G = get_group_elements(group_dims)
    
    # A maximal sum-free set in G = Z_2 x Z_2 x Z_4
    S = [(0, 0, 2), (1, 0, 0), (0, 1, 0)]
    
    print(f"Group G = Z_{group_dims[0]} x Z_{group_dims[1]} x Z_{group_dims[2]}")
    print(f"Size of G is |G| = {len(G)}")
    print(f"Proposed maximal sum-free set S = {S}")
    
    # Verify properties of S
    is_sf = is_sum_free(S, G, group_dims)
    print(f"\nIs S sum-free? {is_sf}")
    
    if not is_sf:
        print("The set S is not sum-free.")
        return

    is_max = is_maximal(S, G, group_dims)
    print(f"Is S maximal by inclusion? {is_max}")

    if not is_max:
        print("The set S is not maximal.")
        return

    # Calculate k(S)
    k_S = find_k_S(S, G, group_dims)
    size_S = len(S)
    size_k_S = len(k_S)
    
    print(f"\nSet S has size |S| = {size_S}")
    print(f"Set k(S) = {{g in G | 2g in S}} is: {k_S}")
    print(f"Set k(S) has size |k(S)| = {size_k_S}")
    
    # Check the condition
    condition_holds = size_k_S > 2 * size_S
    
    print("\nChecking the condition |k(S)| > 2|S|:")
    print(f"{size_k_S} > 2 * {size_S} ?")
    print(f"{size_k_S} > {2 * size_S} ?")
    print(f"Result: {condition_holds}")

    if condition_holds:
        print("\nThe condition is satisfied for this group and set.")
        print(f"The smallest size of such a group is |G| = {len(G)}.")
        
if __name__ == "__main__":
    main()
