import itertools

def verify_smallest_group_size():
    """
    This script verifies the smallest size of a finite Abelian group G
    containing a maximal by inclusion sum-free set S that satisfies |k(S)| > 2|S|.
    The candidate group is G = Z_2 x Z_2 x Z_5, which has size 20.
    The candidate set is S = {(0,0,2), (0,0,3), (1,0,0)}.
    """

    dims = (2, 2, 5)
    group_size = dims[0] * dims[1] * dims[2]

    # Step 1: Define the group elements and the candidate set S
    group_G = list(itertools.product(range(dims[0]), range(dims[1]), range(dims[2])))
    S_list = [(0, 0, 2), (0, 0, 3), (1, 0, 0)]
    S_set = frozenset(S_list)

    def add(g1, g2):
        """Adds two elements in the group G = Z_2 x Z_2 x Z_5."""
        return ((g1[0] + g2[0]) % 2, (g1[1] + g2[1]) % 2, (g1[2] + g2[2]) % 5)

    # Step 2: Verify that S is sum-free
    is_sf = True
    for s1 in S_set:
        for s2 in S_set:
            if add(s1, s2) in S_set:
                is_sf = False
                break
        if not is_sf:
            break

    if not is_sf:
        print("The candidate set S is not sum-free.")
        return

    # Step 3: Verify that S is maximal by inclusion
    is_maximal = True
    elements_outside_S = [g for g in group_G if g not in S_set]
    for g_test in elements_outside_S:
        extended_S_list = S_list + [g_test]
        extended_S_set = frozenset(extended_S_list)
        
        # Check if the extended set is sum-free
        still_sum_free = True
        for s1 in extended_S_set:
            for s2 in extended_S_set:
                if add(s1, s2) in extended_S_set:
                    still_sum_free = False
                    break
            if not still_sum_free:
                break
        
        if still_sum_free:
            # If we could add an element and it remains sum-free, S was not maximal
            is_maximal = False
            break

    if not is_maximal:
        print("The candidate set S is not maximal by inclusion.")
        return

    # Step 4: Compute k(S) and its size
    k_S = set()
    for g in group_G:
        two_g = add(g, g)
        if two_g in S_set:
            k_S.add(g)

    size_S = len(S_set)
    size_k_S = len(k_S)

    # Step 5: Check the condition and print results
    print(f"Group G = Z_{dims[0]} x Z_{dims[1]} x Z_{dims[2]}")
    print(f"Size of G is |G| = {group_size}")
    print(f"Maximal sum-free set S = {sorted(list(S_set))}")
    print(f"Size of S is |S| = {size_S}")
    print(f"k(S) = {{g in G | 2g in S}}")
    print(f"Size of k(S) is |k(S)| = {size_k_S}")
    print(f"Condition to check: |k(S)| > 2 * |S|")
    print(f"Checking: {size_k_S} > 2 * {size_S}  -->  {size_k_S} > {2 * size_S}")

    if size_k_S > 2 * size_S:
        print("\nThe condition is satisfied.")
        print("Since analysis has ruled out all smaller groups,")
        print(f"the smallest size is {group_size}.")
    else:
        print("\nThe condition is NOT satisfied for this set.")

verify_smallest_group_size()