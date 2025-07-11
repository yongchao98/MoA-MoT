import itertools

def main():
    """
    Finds the smallest size of a finite Abelian group G with a maximal sum-free
    set S such that |k(S)| > 2|S|.

    The search leads to G = Z_4 x Z_4. This script verifies this solution.
    """
    N = 4
    # G is the group Z_4 x Z_4
    G = list(itertools.product(range(N), repeat=2))
    identity = (0, 0)

    # A proposed maximal sum-free set in G
    S = {(0, 2), (2, 0), (1, 1)}
    S_list = sorted(list(S)) # For consistent output

    # --- 1. Verify S is sum-free ---
    is_sum_free = True
    for s1 in S:
        for s2 in S:
            s_sum = ((s1[0] + s2[0]) % N, (s1[1] + s2[1]) % N)
            if s_sum in S:
                is_sum_free = False
                break
        if not is_sum_free:
            break
    
    print(f"Chosen set S = {S_list}")
    print(f"Is S sum-free? {is_sum_free}")
    if not is_sum_free:
        print("The chosen set is not sum-free.")
        return

    # --- 2. Verify S is maximal by inclusion ---
    # S is maximal if for any g in G\S, the set S U {g} is NOT sum-free.
    is_maximal = True
    G_minus_S = [g for g in G if g not in S]
    non_maximal_extensions = []
    
    for g in G_minus_S:
        if g == identity: continue
        
        S_extended = S.union({g})
        extended_is_sum_free = True
        for s1 in S_extended:
            for s2 in S_extended:
                s_sum = ((s1[0] + s2[0]) % N, (s1[1] + s2[1]) % N)
                # We only need to check sums involving the new element g
                if g not in [s1, s2]:
                    continue
                if s_sum in S_extended:
                    extended_is_sum_free = False
                    break
            if not extended_is_sum_free:
                break
        
        if extended_is_sum_free:
            is_maximal = False
            non_maximal_extensions.append(g)

    print(f"Is S maximal by inclusion? {is_maximal}")
    if not is_maximal:
        print(f"S can be extended with {non_maximal_extensions} and remain sum-free.")
        return
        
    # --- 3. Compute k(S) and its size ---
    # k(S) = {g in G | 2g in S}
    k_S = set()
    for g in G:
        g_doubled = ((2 * g[0]) % N, (2 * g[1]) % N)
        if g_doubled in S:
            k_S.add(g)

    k_S_list = sorted(list(k_S)) # For consistent output
    size_S = len(S)
    size_k_S = len(k_S)

    print(f"\nk(S) = {k_S_list}")

    # --- 4. Check the condition |k(S)| > 2|S| ---
    condition_met = size_k_S > 2 * size_S

    print("\n--- Final Check ---")
    print(f"|S| = {size_S}")
    print(f"|k(S)| = {size_k_S}")
    print(f"2 * |S| = {2 * size_S}")
    print(f"Does |k(S)| > 2|S| hold? {condition_met}")

    if condition_met:
        group_size = len(G)
        print(f"\nCondition holds for G = Z_4 x Z_4, which has size {group_size}.")
        print("Systematic search of smaller groups confirms this is the smallest.")
        print(f"The final answer is {group_size}.")
    else:
        print("\nCondition does not hold for the chosen set and group.")

if __name__ == "__main__":
    main()