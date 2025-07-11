import itertools

def get_all_subsets(base_set):
    """Generates all subsets of a given set."""
    s = list(base_set)
    return [
        set(subset)
        for i in range(len(s) + 1)
        for subset in itertools.combinations(s, i)
    ]

def is_0_valid(s):
    """Checks if the problem is 0-valid."""
    return 0 in s

def is_1_valid(s):
    """Checks if the problem is 1-valid."""
    return 4 in s

def is_affine(s):
    """Checks if the problem is affine for arity 4."""
    # List of known affine sets for symmetric relations of arity 4
    affine_sets = [
        set(), {0, 1, 2, 3, 4}, {0, 2, 4}, {1, 3},
        {0}, {4}, {0, 4}
    ]
    return s in affine_sets

def is_bijunctive(s):
    """Checks if the problem is bijunctive for arity 4."""
    k = 4
    base_set = set(range(k + 1))
    complement_s = base_set.difference(s)
    # S or its complement must be a subset of {0, k}
    return s.issubset({0, k}) or complement_s.issubset({0, k})

def is_horn(s):
    """Checks if the problem is Horn.
    For non-0-valid cases, only emptyset and {4} are Horn."""
    if 0 in s:
        # This case is already covered by is_0_valid, so it's P-time.
        # We only care about identifying NP-hard cases.
        return True 
    return s == set() or s == {4}


def is_dual_horn(s):
    """Checks if the problem is Dual Horn.
    For non-1-valid cases, only emptyset and {0} are Dual Horn."""
    if 4 in s:
        # This case is already covered by is_1_valid, so it's P-time.
        return True
    return s == set() or s == {0}


def main():
    """
    Main function to find the number of NP-hard sets S.
    """
    base_set = {0, 1, 2, 3, 4}
    all_s = get_all_subsets(base_set)
    
    np_hard_sets = []

    for s in all_s:
        is_p_time = (
            is_0_valid(s) or
            is_1_valid(s) or
            is_affine(s) or
            is_bijunctive(s) or
            is_horn(s) or
            is_dual_horn(s)
        )
        
        if not is_p_time:
            np_hard_sets.append(s)

    print(f"Found {len(np_hard_sets)} NP-hard problems.")
    print("The corresponding sets S are:")
    for s in sorted(list(np_hard_sets), key=lambda x: (len(x), tuple(sorted(list(x))))):
        print(s)

if __name__ == "__main__":
    main()
