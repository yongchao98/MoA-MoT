def solve_magma_problem():
    """
    This function determines for which positive integers n, n-cancellability
    implies mediality in a specific type of magma.

    The condition holds for all positive integers n that are not of the form 2^k - 1
    for some integer k >= 1.
    """

    limit = 100
    print(f"Checking for positive integers n from 1 to {limit}.")
    print("The magma being n-cancellable implies mediality for the following values of n:")

    result_numbers = []

    for n in range(1, limit + 1):
        # A number is of the form 2^k - 1 if and only if n + 1 is a power of 2.
        # A positive number m is a power of 2 if and only if (m & (m - 1)) == 0.
        is_power_of_two_minus_one = ((n + 1) > 0) and (((n + 1) & n) == 0)

        if not is_power_of_two_minus_one:
            result_numbers.append(n)

    # Output the results
    # The problem asks to "output each number in the final equation!".
    # As this is a mathematical problem where the solution is a set of numbers,
    # we will print the list of numbers that satisfy the condition.
    print(", ".join(map(str, result_numbers)))

if __name__ == "__main__":
    solve_magma_problem()
