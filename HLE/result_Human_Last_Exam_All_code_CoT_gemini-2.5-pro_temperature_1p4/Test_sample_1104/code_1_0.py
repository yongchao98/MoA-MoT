def solve_proportionality_problem():
    """
    Solves for the smallest preference profile sizes s1 and s2.

    This is a mathematical problem, and the code prints the derivation and the result.
    """
    k = 100

    print("Step-by-step derivation:")
    print("-------------------------")

    # --- Part 1: Proportional Justified Representation (PJR) ---
    print("\n1. Calculating s1 (for PJR):")
    print(f"   A committee W must leave voter 1 unsatisfied, so W does not contain any of { 'a', 'b', 'c', 'x' }.")
    print(f"   A PJR violation occurs for an unsatisfied group N' if |N'| >= n/k.")
    print(f"   Consider the group N' = {{1}}. This group is unsatisfied and 1-cohesive (|A(1)|=4 >= 1).")
    print(f"   For PJR to be satisfied, this group must NOT cause a violation. So, we need |N'| < n/k.")
    print(f"   Plugging in values: 1 < n / {k}")
    print(f"   This means n > {k}.")
    s1 = k + 1
    print(f"   The smallest integer n greater than {k} is {s1}.")
    print(f"   Therefore, s1 = {s1}.")

    # --- Part 2: Extended Justified Representation (EJR) ---
    print("\n2. Calculating s2 (for EJR):")
    print(f"   An EJR violation occurs for an unsatisfied L-cohesive group N' if |N'| >= L * n/k.")
    print(f"   To satisfy EJR, we must have |N'| < L * n/k for ALL L for which N' is L-cohesive.")
    print(f"   Consider again the group N' = {{1}}. It is L-cohesive for L in {{1, 2, 3, 4}}.")
    print(f"   We must satisfy the condition for all these L values:")
    # L=1
    n_min_l1 = k * 1 / 1
    print(f"   - For L=1: 1 < 1 * n / {k}  => n > {n_min_l1}")
    # L=2
    n_min_l2 = k * 1 / 2
    print(f"   - For L=2: 1 < 2 * n / {k}  => n > {n_min_l2}")
    # L=3
    n_min_l3 = k * 1 / 3
    print(f"   - For L=3: 1 < 3 * n / {k}  => n > {n_min_l3:.2f}")
    # L=4
    n_min_l4 = k * 1 / 4
    print(f"   - For L=4: 1 < 4 * n / {k}  => n > {n_min_l4}")
    print(f"   To satisfy all conditions, n must be greater than the maximum of these lower bounds.")
    s2 = k + 1
    print(f"   The strictest condition is n > {k}. The smallest integer n is {s2}.")
    print(f"   Therefore, s2 = {s2}.")

    # --- Final Answer ---
    print("\n-------------------------")
    print("The final solution (s1, s2) is:")
    print(f"({s1}, {s2})")


solve_proportionality_problem()