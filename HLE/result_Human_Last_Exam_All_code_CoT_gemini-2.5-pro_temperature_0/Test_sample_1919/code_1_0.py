import collections
from itertools import chain, combinations

def solve():
    """
    Determines how many subsets S of {0,1,2,3,4} make the corresponding
    4-variable constraint satisfaction problem NP-hard.
    """
    n = 4
    universe = set(range(n + 1))

    # Helper functions to check for P-time properties based on S
    def is_0_valid(s):
        return 0 in s

    def is_1_valid(s):
        return 4 in s

    def is_affine(s):
        return s in [set(), {0, 2, 4}, {1, 3}, {0, 1, 2, 3, 4}]

    def is_bijunctive(s):
        s_comp = universe - s
        # An interval of length at most 1 is {}, {i}, or {i, i+1}
        intervals = [set()]
        for i in range(n + 1):
            intervals.append({i})
        for i in range(n):
            intervals.append({i, i + 1})
        return s in intervals or s_comp in intervals

    def is_horn(s):
        if not s: return True
        for i in s:
            for j in s:
                lower = max(0, i + j - n)
                upper = min(i, j)
                for k in range(lower, upper + 1):
                    if k not in s:
                        return False
        return True

    def is_dual_horn(s):
        if not s: return True
        for i in s:
            for j in s:
                lower = max(i, j)
                upper = min(n, i + j)
                for k in range(lower, upper + 1):
                    if k not in s:
                        return False
        return True

    # Generate all 2^5 = 32 subsets of {0, 1, 2, 3, 4}
    all_subsets = []
    for r in range(len(universe) + 1):
        for subset in combinations(universe, r):
            all_subsets.append(set(subset))

    p_time_sets = []
    np_hard_sets = []

    for s in all_subsets:
        # Check if S satisfies any of the P-time conditions
        if (is_0_valid(s) or is_1_valid(s) or is_affine(s) or
                is_bijunctive(s) or is_horn(s) or is_dual_horn(s)):
            p_time_sets.append(s)
        else:
            np_hard_sets.append(s)

    print("According to Schaefer's Dichotomy Theorem, the problem is NP-hard if and only if")
    print("the constraint defined by S is not 0-valid, 1-valid, affine, bijunctive, Horn, or dual Horn.")
    print("\nWe test all 2^5 = 32 possible sets S.")
    
    total_sets = len(all_subsets)
    p_time_count = len(p_time_sets)
    np_hard_count = len(np_hard_sets)

    print(f"\nNumber of sets S resulting in a P-time problem: {p_time_count}")
    print(f"Number of sets S resulting in an NP-hard problem: {np_hard_count}")

    print("\nThe calculation is as follows:")
    print(f"Total possible sets S = {total_sets}")
    print(f"Number of P-time sets S = {p_time_count}")
    print(f"Number of NP-hard sets S = {total_sets} - {p_time_count} = {np_hard_count}")

    if np_hard_sets:
        print("\nThe set(s) S for which the problem is NP-hard is/are:")
        for s in np_hard_sets:
            # Print sorted list for consistent output
            print(f"S = {sorted(list(s))}")

solve()