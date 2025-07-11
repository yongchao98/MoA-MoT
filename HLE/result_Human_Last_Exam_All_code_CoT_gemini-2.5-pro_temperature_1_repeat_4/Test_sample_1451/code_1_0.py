def solve_trials():
    """
    Calculates and prints the minimum number of trials T(n) for given values of n.
    """

    def T(n):
        """
        Calculates T(n) based on whether n is even or odd.
        T(n) = 2n - 1 if n is even.
        T(n) = 2n - 2 if n is odd.
        """
        print(f"Calculating for n = {n}:")
        if n % 2 == 0:
            # n is even
            result = 2 * n - 1
            print(f"T({n}) = 2 * {n} - 1 = {result}")
            return result
        else:
            # n is odd
            result = 2 * n - 2
            print(f"T({n}) = 2 * {n} - 2 = {result}")
            return result

    # List of n values to test
    n_values = [2, 3, 1234, 6712]
    results = []

    # Calculate T(n) for each value in the list
    for n_val in n_values:
        results.append(T(n_val))
        print("-" * 20)

    # Print the final comma-separated list of values
    final_answer = ",".join(map(str, results))
    print(f"The values of T(2), T(3), T(1234), and T(6712) are: {final_answer}")
    return final_answer

# Execute the function and capture the final answer for the specified format.
final_answer_string = solve_trials()

# The final answer format required by the prompt
# print(f"<<<{final_answer_string}>>>")