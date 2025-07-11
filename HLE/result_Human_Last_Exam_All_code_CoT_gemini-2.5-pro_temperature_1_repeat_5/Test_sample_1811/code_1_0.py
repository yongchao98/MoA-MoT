import sympy

def solve_task():
    """
    This function provides the formula for the least number of zeros a vector field
    can have on a compact manifold M with a non-empty boundary.
    """

    # Define symbolic variables for the Euler characteristics.
    chi_M = sympy.Symbol('chi(M)')
    chi_dM = sympy.Symbol('chi(partial M)')

    # The formula for the least number of zeros involves these variables.
    # The coefficients/numbers in the formula are 1 and 2.
    c1 = 1
    c2 = 2

    # Construct the formula symbolically.
    # Note: In the actual formula, the division is exact as
    # chi(M) - chi(partial M)/2 is always an integer.
    expression = chi_M - (sympy.Rational(c1, c2)) * chi_dM
    min_zeros_formula = sympy.Abs(expression)

    # Print the explanation and the final formula.
    print("Let chi(M) be the Euler characteristic of the manifold M.")
    print("Let chi(partial M) be the Euler characteristic of its boundary.")
    print("\nThe least number of zeros a vector field can have on M is given by the formula:")
    print(str(min_zeros_formula))

    print("\nThe numbers involved in this final equation are:")
    print(f"The implicit coefficient of chi(M): {c1}")
    print(f"The divisor of chi(partial M): {c2}")

solve_task()