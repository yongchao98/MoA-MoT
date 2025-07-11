import itertools

def is_t_valid(s):
    """Checks if the relation is T-valid (satisfied by all True)."""
    return 4 in s

def is_f_valid(s):
    """Checks if the relation is F-valid (satisfied by all False)."""
    return 0 in s

def is_horn(s):
    """Checks if the relation is Horn-satisfiable."""
    for k in s:
        for l in s:
            lower = max(0, k + l - 4)
            upper = min(k, l)
            for m in range(lower, upper + 1):
                if m not in s:
                    return False
    return True

def is_dual_horn(s):
    """Checks if the relation is dual-Horn-satisfiable."""
    for k in s:
        for l in s:
            lower = max(k, l)
            upper = min(4, k + l)
            for m in range(lower, upper + 1):
                if m not in s:
                    return False
    return True

def is_affine(s):
    """Checks if the relation is affine-satisfiable."""
    affine_sets = [set(), {1, 3}, {0, 2, 4}, {0, 1, 2, 3, 4}]
    return s in affine_sets

def main():
    """
    Main function to iterate through all possible sets S and classify them.
    """
    base_set = {0, 1, 2, 3, 4}
    p_time_sets = []
    np_hard_sets = []
    total_sets = 0

    # Generate all 2^5 = 32 subsets of {0, 1, 2, 3, 4}
    for i in range(len(base_set) + 1):
        for s_tuple in itertools.combinations(base_set, i):
            s = set(s_tuple)
            total_sets += 1

            if (is_t_valid(s) or
                is_f_valid(s) or
                is_horn(s) or
                is_dual_horn(s) or
                is_affine(s)):
                p_time_sets.append(s)
            else:
                np_hard_sets.append(s)

    print("--- Analysis of Constraint Satisfaction Problem Hardness ---")
    print(f"Total number of possible sets S: {total_sets}")
    print(f"Number of sets S for which the problem is in P: {len(p_time_sets)}")
    print(f"Number of sets S for which the problem is NP-hard: {len(np_hard_sets)}")
    
    print("\nThe following 6 sets S result in an NP-hard problem:")
    for s in sorted(list(np_hard_sets), key=lambda x: (len(x), tuple(sorted(list(x))))):
        print(f"S = {s if s else '{}'}")

    print("\nFinal calculation:")
    print(f"Total sets ({total_sets}) - P-time sets ({len(p_time_sets)}) = NP-hard sets ({len(np_hard_sets)})")
    print(f"{total_sets} - {len(p_time_sets)} = {len(np_hard_sets)}")


if __name__ == "__main__":
    main()
