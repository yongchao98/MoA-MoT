import itertools

def is_horn(s):
    """Checks if S is of the form {0, 1, ..., k} or empty."""
    if not s:
        return True
    # It must contain 0 and be a contiguous block from 0.
    return 0 in s and s == set(range(max(s) + 1))

def is_dual_horn(s):
    """Checks if S is of the form {k, k+1, ..., 4} or empty."""
    if not s:
        return True
    # It must contain 4 and be a contiguous block up to 4.
    return 4 in s and s == set(range(min(s), 5))

def is_affine(s):
    """Checks if S is one of the 7 symmetric affine sets."""
    affine_sets = [
        set(), {0}, {4}, {0, 4}, {1, 3}, {0, 2, 4}, {0, 1, 2, 3, 4}
    ]
    return s in affine_sets

def is_bijunctive(s):
    """Checks if S is one of the 7 symmetric bijunctive sets."""
    bijunctive_sets = [
        set(), {0}, {4}, {0, 1}, {3, 4}, {0, 4}, {0, 1, 2, 3, 4}
    ]
    return s in bijunctive_sets

def main():
    """
    Counts the number of NP-hard problems by checking all possible sets S
    against Schaefer's criteria for P-time problems.
    """
    base_set = {0, 1, 2, 3, 4}
    all_s_list = []
    for i in range(len(base_set) + 1):
        for combo in itertools.combinations(base_set, i):
            all_s_list.append(set(combo))

    p_time_count = 0
    np_hard_count = 0
    total_sets = len(all_s_list)
    np_hard_sets = []

    for s in all_s_list:
        # A problem is in P if ANY of the Schaefer conditions are met.
        is_p = (
            0 in s or              # 0-valid
            4 in s or              # 1-valid
            is_horn(s) or
            is_dual_horn(s) or
            is_affine(s) or
            is_bijunctive(s)
        )
        
        if is_p:
            p_time_count += 1
        else:
            np_hard_count += 1
            np_hard_sets.append(s)
    
    print(f"For a 4-variable clause, there are {total_sets} possible satisfaction sets S.")
    print("-" * 20)
    print(f"Number of sets S for which the problem is in P: {p_time_count}")
    print(f"Number of sets S for which the problem is NP-hard: {np_hard_count}")
    print("The NP-hard sets are:", sorted(list(s) for s in np_hard_sets))
    print("-" * 20)
    print("Final Calculation:")
    print(f"{total_sets} (Total Sets) - {p_time_count} (P-time Sets) = {np_hard_count} (NP-hard Sets)")


if __name__ == "__main__":
    main()