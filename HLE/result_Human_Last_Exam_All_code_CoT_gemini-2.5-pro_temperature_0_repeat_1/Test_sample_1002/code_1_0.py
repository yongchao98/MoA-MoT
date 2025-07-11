def solve_limit_problem():
    """
    This function prints the solution to the mathematical problem.
    The problem asks for the limit L = lim_{m->inf} [ln(f(m))/ln(m)] for a given integer k >= 2.

    The derivation shows that this limit is equal to the exponent 'p'
    in the asymptotic relation f(m) ~ C * m^p.

    The value of p was found to be 1 - 1/(2*k).
    """

    # The final formula is L = 1 - 1/(2*k)
    # As requested, we define each number in the final equation.
    numerator = 1
    denominator_const = 2

    # The variable 'k' is given in the problem as an integer k >= 2.
    # We will represent it as a string 'k' in the output formula.
    k_variable_str = 'k'

    print("The derived formula for the limit L is:")
    print(f"L = {numerator} - {numerator} / ({denominator_const} * {k_variable_str})")
    print("\nThis formula is valid for any integer k >= 2.")

    # Example for k=2
    k_val_2 = 2
    limit_val_2 = 1 - 1 / (2 * k_val_2)
    print(f"\nFor k = {k_val_2}, the limit is: 1 - 1/({2*k_val_2}) = {limit_val_2}")

    # Example for k=3
    k_val_3 = 3
    limit_val_3 = 1 - 1 / (2 * k_val_3)
    print(f"For k = {k_val_3}, the limit is: 1 - 1/({2*k_val_3}) = {limit_val_3:.4f}")

solve_limit_problem()