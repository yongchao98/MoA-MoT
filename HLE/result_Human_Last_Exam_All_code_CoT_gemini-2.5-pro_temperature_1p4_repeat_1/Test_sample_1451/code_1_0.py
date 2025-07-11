def solve_and_print():
    """
    Calculates T(n) for a list of n values and prints the detailed equations.
    T(n) is the minimum number of weighings to decide if we have an equal number
    of real (100g) and fake (95g) golden bars among 2n bars.
    """

    def calculate_T_and_print_equation(n):
        """Calculates T(n) based on the derived formula and prints the steps."""
        if n <= 0:
            print(f"T({n}) = 0 (Invalid input)")
            return 0
        
        # Based on the derived logic:
        # T(1) = 1
        # T(n) = 2n - 1, for n even and n >= 2
        # T(n) = 2n - 2, for n odd and n >= 3
        
        if n == 1:
            result = 1
            print(f"T({n}) = {result}")
            return result
        elif n % 2 == 0:  # n is even and >= 2
            result = 2 * n - 1
            print(f"For n={n} (even): T({n}) = 2 * {n} - 1 = {result}")
            return result
        else:  # n is odd and >= 3
            result = 2 * n - 2
            print(f"For n={n} (odd): T({n}) = 2 * {n} - 2 = {result}")
            return result

    values_to_calculate = [2, 3, 1234, 6712]
    final_results = []
    
    print("Calculating T(n) for each value:")
    for n in values_to_calculate:
        res = calculate_T_and_print_equation(n)
        final_results.append(str(res))
    
    print("\nThe comma-separated values for T(2), T(3), T(1234), T(6712) are:")
    print(",".join(final_results))

# Execute the function to see the output
solve_and_print()