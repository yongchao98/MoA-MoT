def solve():
    """
    Calculates T(n), the minimum number of weighings needed for the golden bar problem,
    for a given list of n values.
    The formula derived is T(n) = 2n - 1.
    """
    n_values = [2, 3, 1234, 6712]
    results = []

    def T(n):
        """Calculates the minimum trials T(n) = 2n - 1."""
        return 2 * n - 1

    print("Calculating the results based on the formula T(n) = 2n - 1:\n")
    for n in n_values:
        result = T(n)
        # As requested, output each number in the final equation
        print(f"T({n}) = 2 * {n} - 1 = {result}")
        results.append(result)

    # Format the final answer string as required for the final output block.
    final_answer_string = ",".join(map(str, results))
    print("\n---")
    print("The final values for T(2), T(3), T(1234), and T(6712) are:")
    print(final_answer_string)


solve()