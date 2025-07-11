import sys

def solve():
    """
    This function counts the number of sets S for which the given
    constraint satisfaction problem is NP-hard, based on Schaefer's
    Dichotomy Theorem.
    """

    # These are the sets S that are known to result in P-time problems
    # for symmetric constraints, aside from the 0-valid and 1-valid cases.
    affine_sets = {frozenset(), frozenset({0, 1, 2, 3, 4}), frozenset({0, 2, 4}), frozenset({1, 3})}
    bijunctive_sets = {frozenset(), frozenset({0, 1, 2, 3, 4}), frozenset({0}), frozenset({4}), frozenset({0, 4})}

    # A helper function to check if a set is Horn (an initial segment {0, ..., k})
    def is_horn(s):
        if not s:
            return True
        # To be an initial segment, the set must be equal to {0, 1, ..., max(s)}
        return s == set(range(max(s) + 1))

    # A helper function to check if a set is dual-Horn (a terminal segment {k, ..., 4})
    def is_dual_horn(s):
        if not s:
            return True
        # To be a terminal segment, the set must be equal to {min(s), ..., 4}
        return s == set(range(min(s), 5))

    np_hard_sets = []
    p_time_sets = []

    # Iterate through all 2^5 = 32 possible subsets of {0, 1, 2, 3, 4}
    for i in range(32):
        s = frozenset({j for j in range(5) if (i >> j) & 1})

        # Check for P-time properties
        is_p_time = False
        reason = ""

        # 1. 0-valid
        if 0 in s:
            is_p_time = True
            reason = "0-valid"
        # 2. 1-valid
        elif 4 in s:
            is_p_time = True
            reason = "1-valid"
        # 3. Horn
        elif is_horn(s):
            is_p_time = True
            reason = "Horn"
        # 4. Dual-Horn
        elif is_dual_horn(s):
            is_p_time = True
            reason = "Dual-Horn"
        # 5. Affine
        elif s in affine_sets:
            is_p_time = True
            reason = "Affine"
        # 6. Bijunctive
        elif s in bijunctive_sets:
            is_p_time = True
            reason = "Bijunctive"

        if is_p_time:
            p_time_sets.append(s)
        else:
            np_hard_sets.append(s)

    total_sets = 2**5
    num_p_time = len(p_time_sets)
    num_np_hard = len(np_hard_sets)

    print(f"Out of {total_sets} total possible sets S:")
    print(f"- {num_p_time} sets lead to a problem in P.")
    print(f"- {num_np_hard} sets lead to an NP-hard problem.")
    print("\nThe final equation is:")
    print(f"{total_sets} (Total) - {num_p_time} (P-time) = {num_np_hard} (NP-hard)")


solve()