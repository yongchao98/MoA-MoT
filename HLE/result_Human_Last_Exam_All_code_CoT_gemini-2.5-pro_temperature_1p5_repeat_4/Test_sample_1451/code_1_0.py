def calculate_T_values():
    """
    Calculates T(n) for a given list of n values and prints the steps.
    """
    ns = [2, 3, 1234, 6712]
    results = []

    print("--- Calculation Steps ---")
    for n in ns:
        if n <= 0:
            result = 0
            # No equation needed for non-positive n.
        elif n == 1:
            # For n=1, T(n) is 1.
            result = 1
            print(f"T({n}) = {result}")
        elif n % 2 == 0:
            # For even n >= 2, the formula is T(n) = 2n - 1.
            result = 2 * n - 1
            print(f"T({n}) = 2 * {n} - 1 = {result}")
        else: # n is odd and >= 3
            # For odd n >= 3, the formula is T(n) = 2n - 2.
            result = 2 * n - 2
            print(f"T({n}) = 2 * {n} - 2 = {result}")
        results.append(str(result))

    final_answer = ",".join(results)
    print("\n--- Final Answer ---")
    print("The values of T(2), T(3), T(1234), and T(6712), separated by a comma are:")
    print(final_answer)

calculate_T_values()