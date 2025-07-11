import itertools

def solve_and_demonstrate():
    """
    This script demonstrates the logic for solving the puzzle by checking
    a candidate formula against the given conditions for n=4.
    """
    
    # We choose n=4 as a representative case since n >= 2.
    n = 4
    
    # According to the problem, phi has exactly n distinct atomic variables.
    # We will verify that a simple function depending on only one of these
    # can satisfy the conditions.
    
    # Condition 1: The formula must be true for 2^(n-1) assignments.
    # Here we show all numbers in the equation 2^(n-1) = result
    power = n - 1
    required_true_assignments = 2**power
    
    # Let's consider a formula phi that is logically equivalent to just p1.
    # We will now count how many truth assignments make `p1` true,
    # considering the full set of n=4 variables.
    all_assignments = list(itertools.product([False, True], repeat=n))
    
    # An assignment makes our candidate phi true if the value for p1 is True.
    # In our assignment tuples, p1 is the first element (index 0).
    true_assignment_count = 0
    for assignment in all_assignments:
        if assignment[0] is True:
            true_assignment_count += 1
            
    # --- Verification Step ---
    print(f"For n={n}, the number of required true assignments is given by the equation:")
    print(f"2^({n}-1) = 2^{power} = {required_true_assignments}")
    print(f"\nA formula equivalent to 'p1' results in {true_assignment_count} true assignments.")
    
    if true_assignment_count == required_true_assignments:
        print("This matches the required number, so our candidate formula is valid.")
        
        # --- Conclusion Step ---
        print("\nSince a formula can satisfy the conditions while only depending on one variable (`p1`),")
        print("an equivalent formula `psi` must also depend on `p1`.")
        print("A formula cannot depend on a variable without containing it, so `psi` must have at least one variable.")
        print("The formula `psi = p1` itself uses one variable, showing that 1 is achievable.")
        
        print("\nTherefore, the minimum number of distinct atomic variables required is:")
        print(1)
    else:
        print("Verification failed, indicating an error in the logical model.")

solve_and_demonstrate()