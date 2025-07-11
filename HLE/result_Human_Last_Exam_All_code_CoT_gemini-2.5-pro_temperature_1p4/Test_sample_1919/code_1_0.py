def solve_csp_hardness():
    """
    This function determines the number of NP-hard cases for the given
    constraint satisfaction problem based on Schaefer's Dichotomy Theorem.

    An NP-hard case corresponds to a set S which is a non-empty subset
    of {1, 2, 3} but is not the set {1, 3}.
    """
    
    # The set of potential weight counts that do not contain 0 or 4.
    base_set = {1, 2, 3}
    
    # Total number of subsets of {1, 2, 3} is 2^3.
    num_subsets = 2**len(base_set)
    
    # We exclude two specific P-time cases from these subsets:
    # 1. The empty set: S = {}
    # 2. The affine set: S = {1, 3}
    num_p_time_cases_to_exclude = 2

    # The number of NP-hard sets is the total number of subsets of {1, 2, 3}
    # minus the number of P-time cases within that group.
    num_np_hard_sets = num_subsets - num_p_time_cases_to_exclude
    
    # To satisfy the prompt "output each number in the final equation",
    # we explicitly print the calculation.
    print("The number of NP-hard sets is calculated by considering all subsets of {{1, 2, 3}} and removing the P-time cases.")
    print(f"Total subsets of {{1, 2, 3}}: 2^{len(base_set)} = {num_subsets}")
    print(f"P-time cases to exclude: {num_p_time_cases_to_exclude} (the empty set and the set {{1, 3}})")
    print("\nFinal calculation:")
    print(f"{num_subsets} - {num_p_time_cases_to_exclude} = {num_np_hard_sets}")
    print(f"\nThus, there are {num_np_hard_sets} sets S for which the problem is NP-hard.")

solve_csp_hardness()