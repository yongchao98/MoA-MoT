def solve_weighing_puzzle():
    """
    Calculates and prints the minimum number of trials T(n) for a list of n values.
    The derived formula for T(n) is 2n - 1.
    """
    n_values = [2, 3, 1234, 6712]
    results = []

    print("The minimum number of trials T(n) is given by the formula: T(n) = 2n - 1.")
    print("-" * 30)

    for n in n_values:
        # Calculate T(n) using the formula
        result = 2 * n - 1
        results.append(result)
        # Output each number in the final equation as requested
        print(f"For n = {n}:")
        print(f"T({n}) = 2 * {n} - 1 = {result}")
        print("-" * 30)

    # Format the final answer as a comma-separated string
    final_answer_string = ", ".join(map(str, results))
    print("The values for T(2), T(3), T(1234), and T(6712) are:")
    print(final_answer_string)

# Execute the function
solve_weighing_puzzle()