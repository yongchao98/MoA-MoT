import math

def is_interval(s):
    """Checks if a set of numbers forms a consecutive interval."""
    if not s:
        return True
    s_list = sorted(list(s))
    return s_list[-1] - s_list[0] + 1 == len(s_list)

def main():
    """
    Calculates the number of NP-hard sets S for the given problem.
    """
    base_set = {0, 1, 2, 3, 4}
    num_total_sets = 2**len(base_set)
    
    p_time_sets = set()

    # Category 1: Intervals
    num_intervals = 0
    for i in range(num_total_sets):
        # Generate subset from the bits of i
        s = {elem for j, elem in enumerate(base_set) if (i >> j) & 1}
        if is_interval(s):
            num_intervals += 1
            p_time_sets.add(frozenset(s))

    # Category 2: Affine sets
    affine_sets_tuples = [
        (),          # empty set
        (0,1,2,3,4), # all
        (0,2,4),
        (1,3)
    ]
    affine_sets = {frozenset(t) for t in affine_sets_tuples}
    num_affine = len(affine_sets)
    
    # Union of P-time sets
    for s in affine_sets:
        p_time_sets.add(s)

    num_p_time = len(p_time_sets)
    num_np_hard = num_total_sets - num_p_time

    print(f"Total number of possible sets S is 2^{len(base_set)} = {num_total_sets}.")
    print(f"Number of sets S that are intervals: {num_intervals}")
    print(f"Number of sets S that are affine: {num_affine}")
    print(f"The total number of P-time sets is the size of the union of these two families.")
    # Calculate intersection size for the explanation
    intersection_count = 0
    for s in affine_sets:
        if is_interval(s):
            intersection_count += 1
    
    print(f"|Union| = |Intervals| + |Affine| - |Intersection| = {num_intervals} + {num_affine} - {intersection_count} = {num_p_time}")
    print(f"Number of NP-hard sets = Total sets - P-time sets")
    print(f"The number of NP-hard sets is {num_total_sets} - {num_p_time} = {num_np_hard}")

main()
<<<14>>>