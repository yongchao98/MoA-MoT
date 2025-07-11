def solve_limit_problem():
    """
    Calculates the value of the limit lim_{m->inf} ln(f(m))/ln(m) for a given integer k >= 2.

    The problem is solved analytically, and the result is a function of k.
    This script computes and displays the result for a specific value of k.
    """
    # The problem is defined for an integer k >= 2.
    # The user can change this value to see the result for a different k.
    k = 4

    if not isinstance(k, int) or k < 2:
        print("Error: k must be an integer greater than or equal to 2.")
        return

    # The analytical solution to the limit is k / (k + 1).
    # We will compute the numerator and denominator for the given k.
    numerator = k
    denominator = k + 1
    
    # Calculate the floating-point value of the limit.
    limit_value = numerator / denominator

    # Output the results in a clear, step-by-step format.
    print(f"The problem is to compute the limit of ln(f(m))/ln(m) as m tends to infinity.")
    print(f"The analytical result of this limit is the expression: k / (k + 1)")
    print("-" * 30)
    print(f"For the given value k = {k}:")
    print(f"The final equation for the limit is: {numerator} / ({k} + 1)")
    print(f"The numerator is: {numerator}")
    print(f"The denominator is: {denominator}")
    print(f"The value of the limit is {numerator}/{denominator} = {limit_value}")

if __name__ == '__main__':
    solve_limit_problem()
