import math

def check_horn_property(S_k, k):
    """
    A relation on k variables where satisfaction depends on the number of True variables being in S_k is Horn
    if for any two satisfying assignments t1, t2, their bitwise AND (t1 & t2) is also a satisfying assignment.
    This can be checked by verifying an interval property on the counts of True variables (weights).
    For any weights w1, w2 in S_k, the weight of (t1 & t2) can be any integer in [max(0, w1+w2-k), min(w1, w2)].
    All these possible resulting weights must also be in S_k.
    """
    S_k_set = frozenset(S_k)
    for w1 in S_k_set:
        for w2 in S_k_set:
            lower = max(0, w1 + w2 - k)
            upper = min(w1, w2)
            for m in range(lower, upper + 1):
                if m not in S_k_set:
                    return False
    return True

def check_anti_horn_property(S_k, k):
    """
    Similarly for anti-Horn, the relation is closed under bitwise OR.
    For any weights w1, w2 in S_k, the weight of (t1 | t2) can be any integer in [max(w1, w2), min(k, w1+w2)].
    All these possible resulting weights must also be in S_k.
    """
    S_k_set = frozenset(S_k)
    for w1 in S_k_set:
        for w2 in S_k_set:
            lower = max(w1, w2)
            upper = min(k, w1 + w2)
            for m in range(lower, upper + 1):
                if m not in S_k_set:
                    return False
    return True

def solve():
    """
    Iterates through all 2^5=32 possible sets S and classifies the corresponding
    constraint satisfaction problem as being in P or NP-hard.
    """
    p_sets_count = 0
    total_sets = 2**5
    all_possible_members = {0, 1, 2, 3, 4}

    # The 7 sets S for which the base relation is affine.
    # The family of relations is affine iff the base relation is affine.
    affine_base_sets = [
        frozenset(), frozenset({0, 1, 2, 3, 4}), frozenset({0}), frozenset({4}),
        frozenset({0, 4}), frozenset({0, 2, 4}), frozenset({1, 3})
    ]

    for i in range(total_sets):
        S = frozenset({j for j in range(5) if (i >> j) & 1})
        is_p_time = False

        # 1. 0-valid or 1-valid property
        # The entire family of constraints is 0-valid or 1-valid iff S contains all possible outcomes.
        if S == all_possible_members:
            is_p_time = True

        # 2. Affine property
        if not is_p_time and S in affine_base_sets:
            is_p_time = True

        # 3. Horn property
        if not is_p_time:
            is_horn_family = True
            for k in range(1, 5):  # number of variables
                for c_T in range(0, 5 - k): # number of True constants
                    # The effective satisfaction set for k variables and c_T True constants
                    S_effective = {s - c_T for s in S}
                    S_k_intersect = {v for v in S_effective if 0 <= v <= k}
                    if not check_horn_property(S_k_intersect, k):
                        is_horn_family = False
                        break
                if not is_horn_family:
                    break
            if is_horn_family:
                is_p_time = True

        # 4. Anti-Horn property
        if not is_p_time:
            is_anti_horn_family = True
            for k in range(1, 5):
                for c_T in range(0, 5 - k):
                    S_effective = {s - c_T for s in S}
                    S_k_intersect = {v for v in S_effective if 0 <= v <= k}
                    if not check_anti_horn_property(S_k_intersect, k):
                        is_anti_horn_family = False
                        break
                if not is_anti_horn_family:
                    break
            if is_anti_horn_family:
                is_p_time = True
        
        if is_p_time:
            p_sets_count += 1
            
    np_hard_count = total_sets - p_sets_count
    
    print(f"Total number of possible sets S: {total_sets}")
    print(f"Number of sets S for which the problem is in P: {p_sets_count}")
    print(f"Number of sets S for which the problem is NP-hard: {total_sets} - {p_sets_count} = {np_hard_count}")

solve()