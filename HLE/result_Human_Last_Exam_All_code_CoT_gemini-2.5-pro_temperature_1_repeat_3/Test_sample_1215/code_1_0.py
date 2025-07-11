import itertools

def solve_logic_problem():
    """
    This function solves the logic problem by demonstrating that a formula
    with just one essential variable can satisfy all the given conditions.
    """
    # Let's choose a sample number of total atomic variables, n.
    # The problem states n >= 2. Let's use n = 4 for demonstration.
    n = 4
    print(f"--- Verifying for a system with n = {n} atomic variables ---")

    # Let our candidate formula 'phi' be the simplest possible proposition that
    # might fit the criteria: one that is true if and only if the first
    # atomic variable (p_0) is true.
    # This formula has only one "essential" variable: p_0.
    # The other n-1 variables (p_1, p_2, ...) are "dummy" variables.
    
    # Generate all 2**n possible truth-value assignments for n variables.
    # An assignment is a tuple of booleans, e.g., (True, False, True, False)
    all_assignments = list(itertools.product([True, False], repeat=n))
    total_assignments = len(all_assignments)
    
    print(f"Total possible truth-value assignments: 2^{n} = {total_assignments}")

    # Find all assignments that make our candidate formula 'phi' true.
    # For our 'phi', this means any assignment where the first element is True.
    satisfying_assignments = [v for v in all_assignments if v[0] is True]
    num_satisfying = len(satisfying_assignments)
    
    print(f"Number of assignments that make phi true: {num_satisfying}")
    
    # --- Verification Step ---
    print("\n--- Checking the conditions from the problem statement ---")

    # 1. It has exactly 2^(n-1) discernible truth-value assignments.
    expected_satisfying = 2**(n-1)
    is_condition1_met = (num_satisfying == expected_satisfying)
    print(f"Condition 1: Number of satisfying assignments must be 2^({n}-1).")
    print(f"Calculation: {num_satisfying} == {expected_satisfying}")
    print(f"Condition 1 met: {is_condition1_met}\n")

    # 2. It is not a tautology.
    is_condition2_met = (num_satisfying < total_assignments)
    print("Condition 2: Must not be a tautology (true for all assignments).")
    print(f"Calculation: {num_satisfying} < {total_assignments}")
    print(f"Condition 2 met: {is_condition2_met}\n")

    # 3. For any two distinct truth-value assignments v1 and v2 that make phi true,
    #    there exists at least one atomic variable p such that v1(p) != v2(p).
    # This is true by definition for any set of *distinct* assignments. We can
    # verify it by checking that no two satisfying assignments are identical.
    is_condition3_met = (len(set(satisfying_assignments)) == len(satisfying_assignments))
    print("Condition 3: All satisfying assignments must be distinct.")
    print("Verification: Check if there are any duplicate assignments in the satisfying set.")
    print(f"Condition 3 met: {is_condition3_met}\n")

    # --- Conclusion ---
    if is_condition1_met and is_condition2_met and is_condition3_met:
        print("--- Conclusion ---")
        print("The candidate formula (equivalent to a single proposition 'p_0') satisfies all conditions.")
        print("The truth value of this formula depends on only ONE essential variable.")
        print("Any logically equivalent formula, psi, must also depend on this single variable.")
        print("A formula cannot depend on zero variables, as it would be a tautology or a contradiction, violating condition 1.")
        print("\nTherefore, the minimum number of distinct atomic variables required is 1.")
    else:
        print("The candidate formula did not meet the conditions.")

solve_logic_problem()