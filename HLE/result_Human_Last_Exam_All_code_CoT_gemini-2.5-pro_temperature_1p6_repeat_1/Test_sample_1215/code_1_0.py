import itertools

def solve():
    """
    Solves the propositional logic puzzle by demonstrating with a concrete example.
    """
    # Let's choose a value for n >= 2, as per the problem statement.
    n = 4
    print(f"Let's demonstrate with n = {n} atomic variables (p1, p2, p3, p4).")

    # Let's consider a simple formula phi = p1.
    # Its truth value only depends on the first variable.
    def phi(assignment):
        # assignment is a tuple of booleans, e.g., (True, False, True, False)
        # phi is true if and only if the first variable (p1) is True.
        return assignment[0]

    # 1. Verify that phi satisfies Condition 1:
    # It has exactly 2^(n-1) discernible truth-value assignments.
    
    # Generate all possible 2^n truth assignments for n variables.
    all_assignments = list(itertools.product([False, True], repeat=n))
    
    # Count how many of them make phi true.
    satisfying_assignments = [v for v in all_assignments if phi(v)]
    satisfying_assignment_count = len(satisfying_assignments)
    
    required_count = 2**(n-1)

    print("\nStep 1: Verifying Condition 1 for our example formula phi = p1")
    print(f"Total possible assignments for {n} variables: 2^{n} = {2**n}")
    # The 'equation' from the condition
    print(f"Number of satisfying assignments found: {satisfying_assignment_count}")
    print(f"Required number by condition: 2^(n-1) = 2^({n}-1) = {required_count}")

    if satisfying_assignment_count == required_count:
        print("Conclusion: Condition 1 is met.")
    else:
        print("Conclusion: Condition 1 is NOT met.")

    # 2. Determine the minimum number of variables required for an equivalent formula.
    # This is equal to the number of essential variables of phi.
    # A variable p_i is essential if there exists an assignment v such that
    # changing p_i's value in v changes the value of phi.

    print("\nStep 2: Finding the number of essential variables")
    essential_variables = []
    for i in range(n):
        is_essential = False
        # To check variable p_i+1 (at index i), we iterate through all base assignments
        for base_assignment_list in itertools.product([False, True], repeat=n-1):
            # Create two full assignments that differ only at index i
            temp = list(base_assignment_list)
            
            # v1 has False at index i
            v1_list = temp[:i] + [False] + temp[i:]
            v1 = tuple(v1_list)

            # v2 has True at index i
            v2_list = temp[:i] + [True] + temp[i:]
            v2 = tuple(v2_list)

            # If flipping the variable changes the formula's value, it's essential
            if phi(v1) != phi(v2):
                is_essential = True
                break # Found proof, no need to check other assignments for this variable
        
        if is_essential:
            essential_variables.append(f"p{i+1}")

    print(f"The set of essential variables is: {set(essential_variables)}")
    min_vars_required = len(essential_variables)
    print(f"\nThe minimum number of distinct atomic variables required in any logically equivalent formula psi is the number of essential variables.")
    print(f"Final Answer: {min_vars_required}")

solve()