import itertools

def is_horn(s):
    """Checks if a set S defines a Horn relation for n=4."""
    if not s:
        return True
    s_list = sorted(list(s))
    for i in s_list:
        for j in s_list:
            lower = max(0, i + j - 4)
            upper = min(i, j)
            for k in range(lower, upper + 1):
                if k not in s:
                    return False
    return True

def is_dual_horn(s):
    """Checks if a set S defines a dual-Horn relation for n=4."""
    if not s:
        return True
    s_list = sorted(list(s))
    for i in s_list:
        for j in s_list:
            lower = max(i, j)
            upper = min(4, i + j)
            for k in range(lower, upper + 1):
                if k not in s:
                    return False
    return True

def is_affine(s):
    """Checks if a set S defines an affine relation for n=4."""
    # From characterization: v0=v2, v1=v3, v4=0
    # vk = 1 if k in s, else 0
    v = [1 if i in s else 0 for i in range(5)]
    if v[4] != 0:
        return False
    if v[0] != v[2]:
        return False
    if v[1] != v[3]:
        return False
    return True

def is_bijunctive(s):
    """Checks if a set S defines a bijunctive relation for n=4."""
    # |S intersect {0,1,2,3}| must be even
    # |S intersect {1,2,3,4}| must be even
    count1 = sum(1 for x in s if x in {0, 1, 2, 3})
    count2 = sum(1 for x in s if x in {1, 2, 3, 4})
    return count1 % 2 == 0 and count2 % 2 == 0

def solve():
    """
    Finds the number of sets S for which the CSP is NP-hard.
    """
    base_set = [1, 2, 3]
    np_hard_sets = []
    
    # We only need to check subsets of {1, 2, 3}
    for i in range(len(base_set) + 1):
        for subset in itertools.combinations(base_set, i):
            s = set(subset)
            
            # An NP-hard case must fail ALL P-time conditions.
            # 0-validity (0 in s) and 1-validity (4 in s) are already failed by construction.
            is_p_time = False
            if is_horn(s):
                is_p_time = True
            elif is_dual_horn(s):
                is_p_time = True
            elif is_affine(s):
                is_p_time = True
            elif is_bijunctive(s):
                is_p_time = True
            
            if not is_p_time:
                np_hard_sets.append(s)

    print("The sets S which result in an NP-hard problem are:")
    for s in np_hard_sets:
        # To print empty set correctly
        if not s:
            print("{}")
        else:
            print(s)
            
    print("\nThe total number of such sets is:")
    print(len(np_hard_sets))

solve()