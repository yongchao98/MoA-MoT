import math

def solve_proportionality_problem():
    """
    This script calculates the smallest preference profile sizes (s1, s2)
    for PJR and EJR under the given conditions.
    """

    # --- Problem Parameters ---
    # Committee size
    k = 100
    # The critical group of voters is N_abc = {1, 2, 3, 4, 5, 6}
    group_size = 6
    # The size of the intersection of their ballots {a, b, c}
    intersection_size = 3

    # --- Part 1: Calculation of s1 for PJR ---
    print("--- Calculating s1 (for Proportional Justified Representation) ---")

    # For PJR, a group N' of size |N'| >= n/k that agrees on a candidate 'c'
    # must have 'c' represented in the committee.
    # The group of the first 6 voters agrees on {a,b,c}.
    # If this rule is triggered, one of {a,b,c} must be in the committee,
    # which would satisfy voter 1. To leave voter 1 unsatisfied, this
    # rule must not be triggered for this group.
    # The condition to avoid is: group_size < n / k
    # n > group_size * k
    s1_threshold = group_size * k
    s1 = s1_threshold + 1

    print(f"The PJR axiom applied to the first {group_size} voters requires a candidate from {{'a','b','c'}}")
    print(f"to be in the committee if: {group_size} >= n / {k}")
    print(f"This is equivalent to: n <= {group_size} * {k}, or n <= {s1_threshold}.")
    print("This would satisfy voter 1, contradicting the problem's requirement.")
    print(f"To avoid this, we must have n > {s1_threshold}.")
    print(f"The smallest integer n is {s1_threshold} + 1.")
    print(f"s1 = {group_size} * {k} + 1 = {s1}\n")

    # --- Part 2: Calculation of s2 for EJR ---
    print("--- Calculating s2 (for Extended Justified Representation) ---")

    # For EJR, a group N' of size |N'| >= l*n/k with |intersection| >= l
    # must have at least l members of their union of ballots in the committee.
    # For our group, |intersection| = 3, so we check the condition for l=3.
    # If triggered, |W ∩ {{a,b,c,x,y,z}}| >= 3.
    # But if voter 1 is unsatisfied, W cannot contain {a,b,c,x}, so at most
    # {y,z} can be in W, making the intersection size at most 2.
    # To avoid this contradiction, the EJR rule for l=3 must not be triggered.
    # The condition to avoid is: group_size < l * n / k
    # n > (group_size * k) / l
    l = intersection_size
    s2_threshold = math.floor(group_size * k / l)
    s2 = s2_threshold + 1

    print(f"The EJR axiom applied to the first {group_size} voters with l = {l} requires at least {l} winning candidates")
    print(f"from {{'a','b','c','x','y','z'}} if: {group_size} >= {l} * n / {k}")
    print(f"This is equivalent to: n <= ({group_size} * {k}) / {l}, or n <= {s2_threshold}.")
    print("This contradicts the fact that at most 2 of these candidates (y,z) can be winners.")
    print(f"To avoid this contradiction, we must have n > {s2_threshold}.")
    print(f"The smallest integer n is {s2_threshold} + 1.")
    print(f"s2 = floor({group_size} * {k} / {l}) + 1 = {s2}\n")

    # --- Final Answer ---
    print("--- Final Answer ---")
    result = (s1, s2)
    print(f"The smallest profile size for PJR is s1 = {s1}.")
    print(f"The smallest profile size for EJR is s2 = {s2}.")
    print(f"The solution pair (s1, s2) is: {result}")
    
    print(f"<<<{result}>>>")

solve_proportionality_problem()