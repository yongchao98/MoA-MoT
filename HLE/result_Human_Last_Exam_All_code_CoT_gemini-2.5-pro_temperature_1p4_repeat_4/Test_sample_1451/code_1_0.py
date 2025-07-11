def solve():
    """
    Calculates T(n), the minimum number of trials needed to decide if we have an
    equal number of real and fake golden bars among 2n bars.
    """

    def T(n):
        """
        Calculates the value of T for a given n based on its parity.
        """
        if n == 1:
            return 1
        elif n % 2 == 0:  # n is even
            return 2 * n - 1
        else:  # n is odd and >= 3
            return 2 * n - 2

    # List of n values to calculate T(n) for
    n_values = [2, 3, 1234, 6712]
    results = []

    # Calculate and print the equation for each n
    for n in n_values:
        result = T(n)
        results.append(str(result))
        if n == 1:
            print(f"T({n}) = 1")
        elif n % 2 == 0:
            print(f"T({n}) = 2 * {n} - 1 = {result}")
        else:
            print(f"T({n}) = 2 * {n} - 2 = {result}")
            
    # Print the final answer as a comma-separated string
    final_answer = ",".join(results)
    print("\nFinal comma-separated values:")
    print(final_answer)

solve()