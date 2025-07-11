import itertools

def is_0_valid(s):
    """Checks if the relation is 0-valid (all-false satisfying)."""
    return 0 in s

def is_1_valid(s, n):
    """Checks if the relation is 1-valid (all-true satisfying)."""
    return n in s

def is_bijunctive(s, n):
    """Checks if the relation is bijunctive."""
    # For a symmetric relation, this holds if |S| or |S_complement| is at most 1.
    # The domain is {0, ..., n}, which has size n+1.
    # |S_complement| = (n+1) - |S|.
    # |S_complement| <= 1  is equivalent to |S| >= n.
    # Here n=4, so |S| <= 1 or |S| >= 4.
    return len(s) <= 1 or len(s) >= n

def is_affine(s, n):
    """Checks if the relation is affine over GF(2)."""
    # This check is specific to n=4.
    # The condition is that v_0=v_2=v_4 and v_1=v_3, where v_k = 1 if k in s else 0.
    if n != 4:
        raise ValueError("This affine check is specific for n=4")
    v = [1 if i in s else 0 for i in range(n + 1)]
    return v[0] == v[2] == v[4] and v[1] == v[3]

def is_horn(s, n):
    """Checks if the relation is Horn."""
    # For a symmetric relation, for any i,j in S, the range
    # [max(0, i+j-n), min(i,j)] must be a subset of S.
    if not s: # The empty set corresponds to an always-false relation, which is Horn.
        return True
    for i in s:
        for j in s:
            k_min = max(0, i + j - n)
            k_max = min(i, j)
            for k in range(k_min, k_max + 1):
                if k not in s:
                    return False
    return True

def is_dual_horn(s, n):
    """Checks if the relation is dual Horn."""
    # A relation is dual Horn iff its "negation-transformed" version is Horn.
    # For a symmetric relation S, this corresponds to S' = {n-k | k in S}.
    s_prime = {n - k for k in s}
    return is_horn(s_prime, n)

def solve_csp_hardness_count():
    """
    Counts the number of subsets S of {0,1,2,3,4} for which the
    corresponding constraint satisfaction problem is NP-hard.
    """
    n = 4
    base_set = set(range(n + 1))
    
    # Generate all 2^5 = 32 possible subsets S
    all_s_sets = [set(c) for i in range(len(base_set) + 1) for c in itertools.combinations(base_set, i)]
    
    np_hard_sets = []
    
    for s in all_s_sets:
        # Check if S falls into any of the P-time categories
        is_p_time = (
            is_0_valid(s) or
            is_1_valid(s, n) or
            is_bijunctive(s, n) or
            is_affine(s, n) or
            is_horn(s, n) or
            is_dual_horn(s, n)
        )
        
        # If it's not in any P-time category, it's NP-hard
        if not is_p_time:
            np_hard_sets.append(s)

    print("The sets S for which the problem is NP-hard are:")
    for s in np_hard_sets:
        print(s)
    
    print("\nThe final count of NP-hard sets is the total number of such sets found.")
    final_count = len(np_hard_sets)
    print(f"Final Count = {final_count}")


if __name__ == '__main__':
    solve_csp_hardness_count()