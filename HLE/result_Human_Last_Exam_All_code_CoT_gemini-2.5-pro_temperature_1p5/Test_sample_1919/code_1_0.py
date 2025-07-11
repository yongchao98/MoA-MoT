import math

def is_horn(s, n=4):
    """
    Checks if a set of allowed weights S corresponds to a Horn relation.
    A relation based on the number of true inputs is Horn if for any two
    satisfying assignments u, v, their bitwise AND (u & v) is also satisfying.
    This translates to a condition on the set of allowed weights S:
    for any i, j in S, all possible weights of the intersection of two vectors
    of weights i and j must also be in S. The range of possible weights for
    the intersection is the integer interval [max(0, i + j - n), min(i, j)].
    """
    s_tuple = tuple(s)
    for i in s_tuple:
        for j in s_tuple:
            lower_bound = max(0, i + j - n)
            upper_bound = min(i, j)
            for k in range(lower_bound, upper_bound + 1):
                if k not in s:
                    return False
    return True

def is_dual_horn(s, universe_set):
    """
    Checks if S corresponds to a Dual-Horn relation.
    This is true if the complement of S, with respect to the universe of
    possible weights {0, 1, 2, 3, 4}, is a Horn set.
    """
    complement_s = universe_set - s
    return is_horn(complement_s)

def is_affine(s):
    """
    Checks if S corresponds to an Affine relation.
    For weight-based constraints on 4 variables, the set S must be one of a
    few specific sets for the relation to be affine. This list is derived
    from the property that the set of satisfying assignments must form an
    affine subspace over GF(2).
    """
    affine_sets = [
        frozenset(), frozenset({0, 1, 2, 3, 4}),
        frozenset({0}), frozenset({4}),
        frozenset({0, 4}), frozenset({1, 3}),
        frozenset({0, 2, 4})
    ]
    return s in affine_sets

def solve():
    """
    Calculates the number of NP-hard problems by applying Schaefer's Dichotomy Theorem.
    """
    num_vars = 4
    universe = set(range(num_vars + 1))
    total_subsets = 2**(num_vars + 1)
    
    p_time_sets = set()

    for i in range(total_subsets):
        # Generate the i-th subset S of {0, 1, 2, 3, 4}
        current_s = frozenset({j for j in range(num_vars + 1) if (i >> j) & 1})

        if is_horn(current_s) or is_dual_horn(current_s, universe) or is_affine(current_s):
            p_time_sets.add(current_s)

    num_p_time = len(p_time_sets)
    num_np_hard = total_subsets - num_p_time

    print(f"Total number of possible sets S: 2^5 = {total_subsets}")
    print(f"Number of sets S for which the problem is in P: {num_p_time}")
    print(f"Number of sets S for which the problem is NP-hard: {total_subsets} - {num_p_time} = {num_np_hard}")

solve()