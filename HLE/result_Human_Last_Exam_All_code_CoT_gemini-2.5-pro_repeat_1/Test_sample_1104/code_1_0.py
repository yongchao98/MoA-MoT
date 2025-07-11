def solve_proportionality_problem():
    """
    Calculates s1 and s2 based on the principles of PJR and EJR.
    """
    k = 100  # Committee size
    
    # --- Part 1: Calculation for s1 (PJR) ---
    print("--- Calculating s1 (Proportional Justified Representation) ---")
    
    # For a committee to satisfy PJR while leaving a voter unsatisfied, the group
    # consisting of that single voter must not trigger a PJR violation.
    # A PJR violation for a group N' occurs if all members are unsatisfied AND |N'| >= n/k.
    # For N'={voter 1}, |N'|=1. To avoid violation, we need 1 < n/k, which means n > k.
    s1 = k + 1
    
    print(f"Given committee size k = {k}.")
    print("To leave one voter unsatisfied under PJR, the total number of voters 'n' must be greater than k.")
    print(f"The smallest integer n is k + 1.")
    print(f"s1 = {k} + 1 = {s1}")
    print("-" * 50)
    
    # --- Part 2: Calculation for s2 (EJR) ---
    print("--- Calculating s2 (Extended Justified Representation) ---")

    # The group of the first six voters, N_6 = {A(1), ..., A(6)}, is 3-cohesive
    # because they all approve candidates {a, b, c}.
    group_size = 6
    group_cohesion_l = 3

    # An EJR violation occurs for this group if it is impossible to satisfy its representation claim.
    # This happens if the group triggers the EJR condition: |N_6| >= l * n / k.
    # 6 >= 3 * n / 100  =>  600 >= 3n  =>  200 >= n.
    # If n <= 200, EJR requires |W ∩ {a,b,c,x,y,z}| >= 3.
    # Since voter 1 is unsatisfied, W cannot contain {a,b,c,x}.
    # The requirement becomes |W ∩ {y,z}| >= 3, which is impossible.
    # To avoid this, n must be large enough to fail the trigger condition.
    # We must have n > 200.
    s2_threshold = (group_size * k) / group_cohesion_l
    s2 = int(s2_threshold) + 1

    print(f"The group of the first 6 voters has size |N'| = {group_size}.")
    print(f"This group is l-cohesive for l = {group_cohesion_l}.")
    print("An unavoidable EJR violation occurs if |N'| >= l * n / k.")
    print(f"This is {group_size} >= {group_cohesion_l} * n / {k}, which means n <= {s2_threshold}.")
    print(f"To prevent this violation, n must be greater than {s2_threshold}.")
    print(f"The smallest integer n is floor({s2_threshold}) + 1.")
    print(f"s2 = {int(s2_threshold)} + 1 = {s2}")
    print("-" * 50)
    
    # --- Final Answer ---
    print(f"The final solution is the pair (s1, s2).")
    print(f"({s1}, {s2})")

    # The final answer in the requested format
    # print(f"<<<({s1}, {s2})>>>")

solve_proportionality_problem()