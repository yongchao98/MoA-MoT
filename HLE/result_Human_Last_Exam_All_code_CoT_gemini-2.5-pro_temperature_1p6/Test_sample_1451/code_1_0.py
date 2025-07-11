def solve_weighing_puzzle():
    """
    Calculates T(n), the minimum trials to decide if we have n real and n fake bars.
    The formula derived is T(n) = 2n - 1 for n >= 1.
    """
    values_of_n = [2, 3, 1234, 6712]
    results = []

    print("The minimum number of trials T(n) is given by the formula: T(n) = 2*n - 1.")
    print("-" * 40)

    for n in values_of_n:
        # Calculate each term of the equation
        term1 = 2
        term2 = n
        term3 = 1
        result = term1 * term2 - term3
        
        # Store the final result
        results.append(str(result))
        
        # Output the full equation for each n
        print(f"For n={n}:")
        print(f"T({n}) = {term1} * {term2} - {term3} = {result}")

    print("-" * 40)
    print("The requested values for T(2), T(3), T(1234), and T(6712) are:")
    print(",".join(results))

solve_weighing_puzzle()