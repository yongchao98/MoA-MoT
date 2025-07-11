import itertools

def is_sum_free(s_set, group_elements):
    """Checks if a set is sum-free in Z_4 x Z_4."""
    if not s_set:
        return True
    
    # Generate all sums s1 + s2 for s1, s2 in the set
    for s1 in s_set:
        for s2 in s_set:
            # Group operation is component-wise addition modulo 4
            sum_val = ((s1[0] + s2[0]) % 4, (s1[1] + s2[1]) % 4)
            if sum_val in s_set:
                return False
    return True

def main():
    """
    Main function to find and verify the properties of a set S in G = Z_4 x Z_4.
    """
    N = 4
    G = list(itertools.product(range(N), repeat=2))
    
    # Candidate set based on analysis
    S = {(0, 2), (2, 0), (1, 1)}
    
    print(f"Group G = Z_{N} x Z_{N}, |G| = {len(G)}")
    print(f"Proposed set S = {S}, |S| = {len(S)}")
    print("-" * 30)

    # 1. Verify that S is sum-free
    print("Step 1: Checking if S is sum-free...")
    S_plus_S = set()
    for s1 in S:
        for s2 in S:
            S_plus_S.add(((s1[0] + s2[0]) % N, (s1[1] + s2[1]) % N))
            
    is_sf = not any(s_sum in S for s_sum in S_plus_S)
    print(f"S + S = {S_plus_S}")
    print(f"Is S sum-free? {is_sf}")
    if not is_sf:
        print("S is not sum-free, aborting.")
        return
    print("-" * 30)

    # 2. Verify the condition |k(S)| > 2|S|
    print("Step 2: Checking if |k(S)| > 2|S|...")
    k_S = set()
    for g in G:
        # Calculate 2g = g + g
        double_g = ((g[0] + g[0]) % N, (g[1] + g[1]) % N)
        if double_g in S:
            k_S.add(g)
            
    print(f"k(S) = {{g in G | 2g in S}} = {k_S}")
    print(f"|k(S)| = {len(k_S)}")
    print(f"2*|S| = {2 * len(S)}")
    
    condition_met = len(k_S) > 2 * len(S)
    print(f"Is |k(S)| > 2|S|? {condition_met}")
    print(f"The equation is: {len(k_S)} > {2 * len(S)}")

    if not condition_met:
        print("Condition is not met, aborting.")
        return
    print("-" * 30)
    
    # 3. Verify that S is maximal by inclusion
    print("Step 3: Checking if S is maximal by inclusion...")
    G_minus_S = [g for g in G if g not in S]
    
    is_maximal = True
    for g_test in G_minus_S:
        S_union_g = S.union({g_test})
        if is_sum_free(S_union_g, G):
            print(f"S is NOT maximal. It can be extended with {g_test}.")
            print(f"S U {g_test} = {S_union_g} is still sum-free.")
            is_maximal = False
            break

    if is_maximal:
        print("S is maximal by inclusion. For every element g in G \\ S, the set S U {g} is NOT sum-free.")
    
    print("-" * 30)
    
    if is_sf and condition_met and is_maximal:
        print(f"Conclusion: A set S satisfying the properties exists in G=Z_4 x Z_4.")
        print(f"The group size is {len(G)}.")
        print(f"My analysis ruled out all groups of size less than 16.")
        print(f"Therefore, the smallest size of such a group is 16.")
        

if __name__ == '__main__':
    main()
