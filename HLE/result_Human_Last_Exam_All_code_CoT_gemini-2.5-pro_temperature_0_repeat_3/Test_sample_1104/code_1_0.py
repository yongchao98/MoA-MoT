import math

def solve_proportionality_problem():
    """
    Calculates the smallest preference profile sizes s1 and s2 based on the
    proportionality rules PJR and EJR.
    """
    k = 100  # Committee size
    cohesive_group_size = 6  # Size of the group {1, 2, 3, 4, 5, 6}

    # --- Part 1: Calculation for s1 (PJR) ---
    print("--- Calculating s1 for Proportional Justified Representation (PJR) ---")
    # To leave voter 1 unsatisfied, the common candidates {a,b,c} of the group
    # {1..6} cannot be in the winning committee W.
    # For PJR to hold, this group must not be large enough to guarantee a representative.
    # The condition is: cohesive_group_size < s1 / k
    print(f"The condition for PJR to hold is: |N'| < n / k")
    print(f"Substituting the values: {cohesive_group_size} < s1 / {k}")
    print(f"This implies: s1 > {cohesive_group_size} * {k}")
    
    s1_lower_bound = cohesive_group_size * k
    print(f"So, s1 > {s1_lower_bound}")
    
    # s1 must be the smallest integer greater than the lower bound.
    s1 = s1_lower_bound + 1
    print(f"The smallest integer s1 is {s1}.\n")

    # --- Part 2: Calculation for s2 (EJR) ---
    print("--- Calculating s2 for Extended Justified Representation (EJR) ---")
    # For the group {1..6}, their common candidates {a,b,c} are not in W.
    # So, |C intersect W| = 0.
    # The EJR condition is: |C intersect W| >= floor(|N'| * k / n)
    # This becomes: 0 >= floor(cohesive_group_size * k / s2)
    print(f"The condition for EJR to hold is: 0 >= floor(|N'| * k / n)")
    print(f"Substituting the values: 0 >= floor({cohesive_group_size} * {k} / s2)")
    
    # For floor(x) to be 0, we need x < 1.
    print(f"This implies: ({cohesive_group_size} * {k}) / s2 < 1")
    s2_lower_bound = cohesive_group_size * k
    print(f"So, {s2_lower_bound} < s2")
    
    # s2 must be the smallest integer greater than the lower bound.
    s2 = s2_lower_bound + 1
    print(f"The smallest integer s2 is {s2}.\n")

    # --- Final Result ---
    print(f"The solution is the pair (s1, s2) = ({s1}, {s2})")
    return (s1, s2)

# Execute the function to print the solution
solve_proportionality_problem()