import math

def solve_proportionality_problem():
    """
    Calculates s1 and s2 based on the definitions of PJR and EJR
    and the given problem parameters.
    """
    # Given parameters from the problem description
    k = 100  # Committee size

    # From the ballots, we identify a key cohesive group.
    # Voters {1,2,3,4,5,6} all approve of candidates {a,b,c}.
    # Let N be this group and C_full be this set of candidates.
    N_size = 6
    
    # The condition that voter 1 is unsatisfied (A(1)={a,b,c,x}) implies
    # that {a,b,c} must be disjoint from the winning committee W.

    print("--- Calculating s_1 (for Proportional Justified Representation) ---")
    
    # PJR Violation Condition: |N| >= n/k
    # To satisfy PJR, we must have: |N| < s_1 / k
    # This prevents the group N from constituting a violation.
    print(f"To satisfy PJR, the size of the profile s_1 must satisfy the inequality:")
    print(f"{N_size} < s_1 / {k}")
    
    # Solving for s_1:
    # s_1 > N_size * k
    s1_lower_bound = N_size * k
    print(f"s_1 > {N_size} * {k}")
    print(f"s_1 > {s1_lower_bound}")
    
    # The smallest integer s_1 greater than the lower bound
    s1 = s1_lower_bound + 1
    print(f"The smallest integer value for s_1 is {s1}.")
    print("\n" + "="*50 + "\n")

    print("--- Calculating s_2 (for Extended Justified Representation) ---")

    # EJR Violation Condition: for some l, |N| >= l * n/k AND |C| = l AND |C intersect W| < l
    # For group N, any non-empty subset C_sub of {a,b,c} has |C_sub intersect W| = 0.
    # So if |C_sub|=l, then 0 < l, so |C_sub intersect W| < l is always true.
    # To satisfy EJR, we must ensure |N| < l * s_2 / k for all these cases.
    
    # The set C={a,b,c} allows for subsets of size l=1, l=2, and l=3.
    possible_l = [1, 2, 3]
    
    print("To satisfy EJR, s_2 must satisfy |N| < l * s_2 / k for all relevant l.")
    print("This is equivalent to s_2 > |N| * k / l.")
    
    constraints = []
    for l in possible_l:
        s2_lower_bound_for_l = (N_size * k) / l
        constraints.append(s2_lower_bound_for_l)
        print(f"\nFor l = {l}:")
        print(f"The condition is {N_size} < {l} * s_2 / {k}")
        print(f"s_2 > {N_size} * {k} / {l}")
        print(f"s_2 > {s2_lower_bound_for_l}")

    # To satisfy all EJR conditions, s_2 must be greater than the maximum of these bounds.
    s2_lower_bound = max(constraints)
    print(f"\nTo satisfy all these conditions simultaneously, s_2 must be greater than the maximum lower bound.")
    print(f"The strictest constraint is s_2 > {s2_lower_bound}.")

    # The smallest integer s_2 greater than this bound
    s2 = math.floor(s2_lower_bound) + 1
    print(f"The smallest integer value for s_2 is {s2}.")
    print("\n" + "="*50 + "\n")

    print(f"The final solution pair is (s_1, s_2).")
    print(f"Result: ({s1}, {s2})")

solve_proportionality_problem()