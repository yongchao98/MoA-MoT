def solve_and_print():
    """
    This script calculates the minimum number of weighings T(n) needed to decide
    if there is an equal number of real and fake golden bars among 2n bars.
    It prints the formula and result for each requested value of n, and then
    prints a comma-separated list of the results.
    """

    def get_t_value_and_equation(n):
        """
        Applies the derived formula to calculate T(n) and returns the value
        along with a string explaining the calculation.
        """
        if not isinstance(n, int) or n <= 0:
            return None, "Error: n must be a positive integer."
        
        # Formula for T(n)
        if n == 1:
            val = 1
            eq = f"For n=1, T(1) = 1"
        elif n % 2 == 0:  # n is even
            val = 2 * n - 1
            eq = f"For even n={n}, T({n}) = 2 * {n} - 1 = {val}"
        else:  # n is odd and > 1
            val = 2 * n - 2
            eq = f"For odd n={n}, T({n}) = 2 * {n} - 2 = {val}"
        return val, eq

    # The values of n for which T(n) is requested
    ns_to_calculate = [2, 3, 1234, 6712]
    
    final_results = []
    
    print("Calculating the minimum number of trials T(n):")
    print("-" * 45)
    
    for n in ns_to_calculate:
        value, equation_string = get_t_value_and_equation(n)
        if value is not None:
            final_results.append(str(value))
            print(equation_string)
            
    print("-" * 45)
    print("The final values for T(2), T(3), T(1234), and T(6712) are:")
    print(",".join(final_results))

# Execute the function to see the output
solve_and_print()