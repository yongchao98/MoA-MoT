import itertools

def solve():
    """
    This function determines how many subsets S of {0,1,2,3,4} lead to an NP-hard problem.
    """
    base_set = {0, 1, 2, 3, 4}
    np_hard_sets = []

    # Generate all 2^5 = 32 possible subsets of base_set
    for i in range(len(base_set) + 1):
        for s_tuple in itertools.combinations(base_set, i):
            s = set(s_tuple)

            # A problem is in P if it is 0-valid, 1-valid, or affine.
            # Otherwise, it is NP-hard.

            # 1. Check for 0-validity (is 0 in S?)
            is_0_valid = 0 in s

            # 2. Check for 1-validity (is 4 in S?)
            is_1_valid = 4 in s

            # 3. Check for affinity (are all numbers in S of the same parity?)
            # An empty set is considered affine.
            if not s:
                is_affine = True
            else:
                has_even = any(x % 2 == 0 for x in s)
                has_odd = any(x % 2 != 0 for x in s)
                is_affine = not (has_even and has_odd)

            # If none of the P-time conditions are met, the problem is NP-hard.
            if not is_0_valid and not is_1_valid and not is_affine:
                np_hard_sets.append(s)

    print("The sets S for which the problem is NP-hard are:")
    for s in np_hard_sets:
        # Print the elements of the set, sorted for consistent output
        print(sorted(list(s)))
    
    print("\nThis gives a total of:")
    print(len(np_hard_sets))

solve()