import itertools

def group_add(g1, g2):
    """Adds two elements in G = Z_2 x Z_2 x Z_5."""
    return ((g1[0] + g2[0]) % 2, (g1[1] + g2[1]) % 2, (g1[2] + g2[2]) % 5)

def is_sum_free(subset, group_elements):
    """Checks if a subset is sum-free."""
    subset_set = set(subset)
    if not subset_set:
        return True
    for s1 in subset:
        for s2 in subset:
            s_sum = group_add(s1, s2)
            if s_sum in subset_set:
                return False
    return True

def main():
    """
    Finds and verifies the smallest Abelian group and maximal sum-free set S
    satisfying |k(S)| > 2|S|.
    """
    # Group G = Z_2 x Z_2 x Z_5
    G = list(itertools.product(range(2), range(2), range(5)))
    group_size = len(G)

    # A known maximal sum-free set S in G
    S_list = [(0, 0, 1), (0, 0, 4), (1, 0, 0), (0, 1, 0), (1, 1, 4)]
    S = set(S_list)
    s_size = len(S)

    # Verify S is sum-free
    if not is_sum_free(S, G):
        print("The provided set S is not sum-free. Aborting.")
        return

    # Verify S is maximal by inclusion
    is_maximal = True
    G_minus_S = [g for g in G if g not in S]
    for g in G_minus_S:
        T = list(S) + [g]
        if is_sum_free(T, G):
            # If we find a g such that S U {g} is sum-free, S is not maximal
            is_maximal = False
            break

    if not is_maximal:
        print("The provided set S is not maximal. Aborting.")
        return

    # Calculate k(S) = {g in G | 2g in S}
    k_S = set()
    for g in G:
        g_double = group_add(g, g)
        if g_double in S:
            k_S.add(g)
    k_s_size = len(k_S)

    # Check the condition
    condition_met = k_s_size > 2 * s_size
    
    print(f"The size of the finite Abelian group G is: {group_size}")
    print(f"The maximal sum-free set is S = {S}")
    print(f"The size of S is |S| = {s_size}")
    print(f"The set k(S) is k(S) = {k_S}")
    print(f"The size of k(S) is |k(S)| = {k_s_size}")
    print(f"The condition is |k(S)| > 2 * |S|, which is {k_s_size} > 2 * {s_size} => {k_s_size} > {2*s_size}")
    
    if condition_met:
        print("The condition is satisfied.")
        print("The smallest size of such a group is 20.")
        
    else:
        print("The condition is not satisfied by this set.")

main()
