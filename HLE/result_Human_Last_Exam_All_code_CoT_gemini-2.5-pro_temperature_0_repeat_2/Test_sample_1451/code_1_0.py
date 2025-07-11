def solve_weighing_puzzle():
    """
    This script calculates the minimum number of trials T(n) needed to decide
    if we have an equal number of real and fake golden bars among 2n bars.
    """

    def get_t_n_with_equation(n):
        """
        Calculates T(n) based on the derived formula and returns the value
        and a string representation of the calculation.
        """
        if n <= 0:
            return 0, f"T({n}) is not defined for n <= 0"
        
        # For n=1, we have 2 bars. One weighing (bar1 vs bar2) is sufficient.
        if n == 1:
            result = 1
            equation = f"T({n}) = {result}"
            return result, equation

        # For n > 1 and n is even, the worst-case formula is 2n - 1.
        if n % 2 == 0:
            result = 2 * n - 1
            equation = f"T({n}) = 2 * {n} - 1 = {result}"
            return result, equation
        # For n > 1 and n is odd, the worst-case formula is 2n - 2.
        else:
            result = 2 * n - 2
            equation = f"T({n}) = 2 * {n} - 2 = {result}"
            return result, equation

    # The specific values of n we need to solve for.
    n_values_to_solve = [2, 3, 1234, 6712]

    final_values = []
    
    print("The minimum number of trials T(n) is calculated as follows:")
    
    for n_val in n_values_to_solve:
        value, eq_str = get_t_n_with_equation(n_val)
        final_values.append(str(value))
        # As requested, outputting each number in the final equation.
        print(eq_str)

    # The final answer, as a comma-separated string of values.
    print("\nThe required values for T(2), T(3), T(1234), and T(6712) are:")
    print(",".join(final_values))

solve_weighing_puzzle()