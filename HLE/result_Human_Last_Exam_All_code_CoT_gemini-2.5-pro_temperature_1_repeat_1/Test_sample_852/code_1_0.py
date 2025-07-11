import itertools

def is_sum_free(s, n1=4, n2=4):
    """Checks if a set s is sum-free in Z_n1 x Z_n2."""
    for e1 in s:
        for e2 in s:
            sum_val = ((e1[0] + e2[0]) % n1, (e1[1] + e2[1]) % n2)
            if sum_val in s:
                return False
    return True

def main():
    """
    Finds the smallest size of a finite Abelian group G containing a maximal
    by inclusion sum-free set S that satisfies |k(S)| > 2|S|.
    The search suggests G = Z_4 x Z_4 is a candidate.
    """
    n1, n2 = 4, 4
    group_size = n1 * n2
    
    # Candidate set S in Z_4 x Z_4
    S = frozenset([(0, 2), (2, 2), (1, 0)])

    # 1. Check if S itself is sum-free.
    if not is_sum_free(S, n1, n2):
        print(f"The set S={S} is not sum-free.")
        return

    # 2. Check if S is maximal by inclusion.
    # Iterate through all elements g not in S.
    all_elements = frozenset(itertools.product(range(n1), range(n2)))
    elements_to_check = all_elements - S - frozenset([(0,0)])

    is_maximal = True
    for g in elements_to_check:
        s_prime = S.union({g})
        if is_sum_free(s_prime, n1, n2):
            is_maximal = False
            print(f"S is not maximal, because S union {{ {g} }} is still sum-free.")
            break
            
    if is_maximal:
        print(f"S = {set(S)} is a maximal by inclusion sum-free set in Z_{n1} x Z_{n2}.")
        
        # 3. Calculate k(S) and check the condition |k(S)| > 2|S|
        K = frozenset(g for g in all_elements if (2*g[0]) % n1 == 0 and (2*g[1]) % n2 == 0)
        Im_phi = frozenset(((2*g[0]) % n1, (2*g[1]) % n2) for g in all_elements)
        
        S_e = S.intersection(Im_phi)
        
        k_S_size = len(S_e) * len(K)
        S_size = len(S)
        
        print(f"|S| = {S_size}")
        print(f"|k(S)| = |S_e|*|K| = {len(S_e)} * {len(K)} = {k_S_size}")
        
        if k_S_size > 2 * S_size:
            print(f"The condition |k(S)| > 2|S| is met ({k_S_size} > 2 * {S_size}).")
            print(f"The smallest size of such a group is {group_size}.")
            final_equation = f"{k_S_size} > 2 * {S_size}"
            print(final_equation)
        else:
            print(f"The condition |k(S)| > 2|S| is NOT met ({k_S_size} <= 2 * {S_size}).")
    else:
        print("The chosen S was not maximal. A different set or group is needed.")

main()
