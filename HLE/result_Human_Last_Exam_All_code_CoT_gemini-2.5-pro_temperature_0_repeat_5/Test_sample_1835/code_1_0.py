def demonstrate_generality_constraint():
    """
    This script models the Generality Constraint by showing how understanding
    a predicate for one object allows understanding it for all objects.
    """
    # Let's define a domain of objects.
    domain = [10, 20, 30, 41, 50]

    # Let 'a' be a specific object we understand a proposition about.
    a = 20

    # Let F(x) be a predicate, for example, "x is a multiple of 10".
    # Understanding F(a) means we have a general grasp of this predicate.
    def F(x):
        return x % 10 == 0

    print("--- Step 1: Understanding the proposition F(a) ---")
    print(f"Let the predicate F(x) be 'x is a multiple of 10'.")
    print(f"Let the object 'a' be the number {a}.")
    result_fa = F(a)
    print(f"The proposition F(a) is F({a}), which evaluates to: {result_fa}\n")

    print("--- Step 2: Combining F(x) with Universal Quantification (∀x) ---")
    print("Assuming we also understand 'for all x' (∀x), we can form the new proposition '∀x F(x)'.")
    print(f"This means checking if F(x) is true for every x in our domain: {domain}.")
    print("The proposition '∀x F(x)' expands to the following equation:")

    # Build and print the equation string with object names
    equation_str = " and ".join([f"F({x})" for x in domain])
    print(equation_str)

    # Build and print the equation string with the result of each F(x)
    result_values = [F(x) for x in domain]
    result_str = " and ".join([str(val) for val in result_values])
    print(f"= {result_str}")

    # Calculate and print the final boolean result
    final_result = all(result_values)
    print(f"= {final_result}")

demonstrate_generality_constraint()