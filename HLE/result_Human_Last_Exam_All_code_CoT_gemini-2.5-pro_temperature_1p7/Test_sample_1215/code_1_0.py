import itertools

def solve():
    """
    Solves the logic puzzle by demonstrating that a formula dependent on a single
    variable can satisfy the given conditions.
    """
    # The problem states n >= 2. Let's choose a sample value for demonstration.
    n = 3
    
    print(f"Let's analyze the problem for a sample case where n = {n}.")
    print(f"The total number of truth-value assignments for {n} variables is 2^{n} = {2**n}.")

    # Condition 1: The formula must have exactly 2^(n-1) models.
    required_models = 2**(n - 1)
    print(f"Condition 1: The number of true assignments (models) must be 2^({n}-1) = {required_models}.")

    # Let's consider a simple formula 'phi' that is logically equivalent to p1.
    # We will check if this formula satisfies the conditions.
    print("\nLet's test a candidate formula 'phi' that depends only on the first variable, p1.")
    
    # Generate all 2^n possible truth assignments for (p1, p2, ..., pn).
    assignments = list(itertools.product([True, False], repeat=n))
    
    # Count the models for our candidate formula (phi is true iff p1 is true).
    models = [v for v in assignments if v[0] is True]
    num_models_found = len(models)

    print(f"For a formula equivalent to p1, the models are the assignments where p1 is True:")
    # for model in models:
    #     print(f"  {model}")
    print(f"The count of such models is {num_models_found}.")

    # Verify if our candidate formula works.
    if num_models_found == required_models:
        print("This matches the required number of models. Condition 1 is satisfied.")
    else:
        print("This does not match the required number. Condition 1 is not satisfied.")

    # Condition 2: Not a tautology.
    # A tautology is true for all 2^n assignments.
    print(f"Condition 2: The formula must not be a tautology. Our formula is true for {num_models_found} out of {2**n} cases, so it is not a tautology.")

    # Condition 3 is always true by definition of distinct assignments.
    print("Condition 3 is satisfied for any formula.")

    print("\nConclusion:")
    print("We have shown that a formula depending on only one variable can satisfy all the problem's conditions.")
    
    # A formula `psi` equivalent to our `phi` would be just `p1`.
    # This requires only one atomic variable.
    min_vars = 1
    
    # The number of variables cannot be 0, as that would be a constant (True or False for all assignments),
    # which does not have 2^(n-1) models for n>=2.
    print(f"An equivalent formula `psi` only needs to depend on that single variable, p1.")
    print(f"Therefore, the minimum number of distinct atomic variables required is {min_vars}.")
    final_answer = 1
    print(f"\nFinal Equation: Minimum variables = {final_answer}")

solve()
<<<1>>>