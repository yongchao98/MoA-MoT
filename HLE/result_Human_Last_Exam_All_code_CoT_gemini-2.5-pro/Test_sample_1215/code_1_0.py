import itertools

def solve_logic_problem():
    """
    Solves the propositional logic problem by deduction and demonstrates the result with code.
    """

    # --- Step 1: Logical Deduction (as explained in the thinking steps) ---

    # The problem asks for the minimum number of variables a formula `psi` can have,
    # where `psi` is equivalent to a formula `phi` that has n variables and 2^(n-1) satisfying assignments.
    # This minimum number is the number of variables `k` that the underlying truth function semantically depends on.

    # If a function on n variables depends on only k variables, its total number of satisfying assignments is:
    # (Number of satisfying assignments on k variables, N_k) * 2^(n-k)
    # We are given this total is 2^(n-1).
    # So, N_k * 2^(n-k) = 2^(n-1)
    # This simplifies to N_k = 2^(k-1).

    # We need to find the minimum k >= 0 for which a function on k variables can have 2^(k-1) models.
    # - For k=0, N_0 = 2^(-1) = 0.5. Impossible for a function to have 0.5 models.
    # - For k=1, N_1 = 2^(0) = 1. A function f(p) = p has exactly 1 model (when p is True). This is possible.

    # Therefore, the minimum number of variables the function must depend on is 1.

    # --- Step 2: Code Demonstration ---
    # We will demonstrate this for a sample n. Let's choose n=4.
    # A formula `phi` would have 4 variables and 2^(4-1) = 8 satisfying assignments.
    # We hypothesize it can depend on k=1 variable.

    n = 4
    k = 1

    print(f"Demonstration for n = {n}")
    print("-------------------------")
    print("A formula phi has n variables and is true for 2^(n-1) assignments.")
    
    required_models = 2**(n - 1)
    print(f"For n={n}, the required number of satisfying assignments is 2^({n}-1) = {required_models}.")

    print(f"\nLet's assume the formula is equivalent to a simpler formula that only depends on k={k} variable(s).")
    print("We test the simplest k=1 function: a function that is true if and only if the first variable is True.")

    # In our k=1 test case (f(p1) = p1), the number of models in the k-space is N_k = 1.
    N_k = 1
    
    # Let's verify our derived equation: N_k * 2^(n-k) = 2^(n-1)
    lhs_val = N_k * (2**(n - k))
    rhs_val = 2**(n - 1)
    
    print("\nVerifying the core equation: N_k * 2^(n-k) = 2^(n-1)")
    print("Substituting the values for our demonstration case (n=4, k=1, N_1=1):")
    print(f"Equation: {N_k} * 2^({n}-{k}) = 2^({n}-1)")
    print(f"Calculation: {N_k} * {2**(n-k)} = {2**(n-1)}")
    print(f"Result: {lhs_val} = {rhs_val}")
    print("The equation holds true.")


    # Now, let's confirm by iterating through all 2^n possibilities.
    satisfying_assignments_count = 0
    variables = [f'p{i+1}' for i in range(n)]
    
    # Generate all 2^n possible truth assignments
    for assignment in itertools.product([False, True], repeat=n):
        # Our simple function only depends on the first variable
        if assignment[0] is True:
            satisfying_assignments_count += 1
            
    print(f"\nSimulation result: Counting satisfying assignments for a function of {n} variables that only depends on the first variable...")
    print(f"The simulation found {satisfying_assignments_count} satisfying assignments.")

    if satisfying_assignments_count == required_models:
        print("This matches the required number. This demonstrates that a formula satisfying the conditions can indeed depend on just 1 variable.")
    else:
        print("This does not match the required number. There is an error in the logic.")
    
    print("\nSince k=1 is possible and k=0 is not, the minimum number of variables required is 1.")

solve_logic_problem()
print("\n<<<1>>>")