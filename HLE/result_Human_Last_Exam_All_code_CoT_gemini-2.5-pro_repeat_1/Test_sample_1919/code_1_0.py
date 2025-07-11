def solve():
    """
    This function counts the number of sets S for which the corresponding
    constraint satisfaction problem is NP-hard.
    
    A set S leads to a problem in P if it falls into one of Schaefer's tractable classes.
    For symmetric constraints on 4 variables, this means S must be:
    - An interval: e.g., {1, 2, 3}.
    - An affine set: one of {}, {1, 3}, {0, 2, 4}, {0, 1, 2, 3, 4}.
    
    Any set S that does not meet either of these criteria results in an NP-hard problem.
    """
    
    affine_sets = {
        frozenset(),
        frozenset({1, 3}),
        frozenset({0, 2, 4}),
        frozenset({0, 1, 2, 3, 4})
    }

    def is_interval(s_set):
        if not s_set:
            return True
        min_val = min(s_set)
        max_val = max(s_set)
        return len(s_set) == (max_val - min_val + 1)

    np_hard_count = 0
    np_hard_sets = []

    # There are 2^5 = 32 possible subsets of {0, 1, 2, 3, 4}
    for i in range(32):
        current_set = set()
        if (i >> 0) & 1: current_set.add(0)
        if (i >> 1) & 1: current_set.add(1)
        if (i >> 2) & 1: current_set.add(2)
        if (i >> 3) & 1: current_set.add(3)
        if (i >> 4) & 1: current_set.add(4)
        
        # Check if the problem is in P
        is_p_time = False
        if is_interval(current_set):
            is_p_time = True
        elif frozenset(current_set) in affine_sets:
            is_p_time = True

        if not is_p_time:
            np_hard_count += 1
            np_hard_sets.append(sorted(list(current_set)))
            
    print(f"There are {np_hard_count} sets S for which the problem is NP-hard.")
    print("These sets are:")
    for s in np_hard_sets:
        print(s)

solve()
<<<14>>>