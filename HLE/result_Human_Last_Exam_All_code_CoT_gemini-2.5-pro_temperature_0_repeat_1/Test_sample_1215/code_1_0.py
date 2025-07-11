import itertools

def solve():
    """
    This script verifies that a formula with 2 essential variables can satisfy
    the conditions of the problem, confirming that 2 is a possible value.
    Our logical derivation shows it is also the minimum possible value.

    Let's test with n=2 and the formula phi = p_1 XOR p_2.
    The variables are p_1 and p_2. Both are essential.
    So, k = n = 2.
    """
    n = 2
    # The set of essential variables, indexed from 0. Here, {0, 1} for p_1, p_2.
    essential_vars_indices = {0, 1}
    k = len(essential_vars_indices)

    print(f"--- Verifying for n = {n} and k = {k} ---")
    print(f"Formula (phi): p_1 XOR p_2")

    # A truth assignment is a tuple, e.g., (True, False) for (p_1, p_2)
    # The formula p_1 XOR p_2 is true if p_1 != p_2
    def evaluate_phi(assignment):
        return assignment[0] != assignment[1]

    # Generate all 2^n possible truth assignments
    all_assignments = list(itertools.product([False, True], repeat=n))

    # Find all assignments that make phi true
    satisfying_assignments = [v for v in all_assignments if evaluate_phi(v)]

    print(f"\nTotal assignments: {2**n}")
    print(f"Satisfying assignments (V_true): {satisfying_assignments}")

    # --- Check Condition 1 ---
    # It has exactly 2^(n-1) discernible truth-value assignments.
    num_satisfying = len(satisfying_assignments)
    required_num_satisfying = 2**(n - 1)
    cond1_met = (num_satisfying == required_num_satisfying)
    print("\n--- Condition 1 Check ---")
    print(f"Number of satisfying assignments: {num_satisfying}")
    print(f"Required number (2^(n-1)): {required_num_satisfying}")
    print(f"Condition 1 Met: {cond1_met}")

    # --- Check Condition 2 ---
    # It is not a tautology.
    cond2_met = (num_satisfying < 2**n)
    print("\n--- Condition 2 Check ---")
    print(f"Is it a tautology? {not cond2_met}")
    print(f"Condition 2 Met: {cond2_met}")

    # --- Check Condition 3 (with the non-trivial interpretation) ---
    # For any two distinct v1, v2 in V_true, they must differ on an essential variable.
    cond3_met = True
    if len(satisfying_assignments) >= 2:
        # Generate all pairs of distinct satisfying assignments
        for v1, v2 in itertools.combinations(satisfying_assignments, 2):
            # Check if they differ on at least one essential variable
            differ_on_essential = False
            for i in essential_vars_indices:
                if v1[i] != v2[i]:
                    differ_on_essential = True
                    break
            if not differ_on_essential:
                cond3_met = False
                break
    print("\n--- Condition 3 Check (Essential Variable Interpretation) ---")
    print(f"Condition 3 Met: {cond3_met}")

    if cond1_met and cond2_met and cond3_met:
        print("\nConclusion: A formula with n=2 variables can satisfy all conditions.")
        print("The number of essential variables in this formula is k=2.")
        print("Our logical argument shows that k must be >= 2.")
        print("\nTherefore, the minimum number of distinct atomic variables required is 2.")
    else:
        print("\nVerification failed. The example does not meet the conditions.")

solve()
<<<2>>>