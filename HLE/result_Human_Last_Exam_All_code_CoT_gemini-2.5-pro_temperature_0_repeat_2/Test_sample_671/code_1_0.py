import sympy

def generate_and_print_polynomial(n):
    """
    Calculates the polynomial f_n(p) using the derived formula and prints the result.

    The formula is f_n(p) = (p^n - (1-p)^n) / (2p - 1).

    Args:
        n (int): The index of the polynomial in the sequence.
    """
    if not isinstance(n, int) or n < 1:
        print("Error: n must be a positive integer.")
        return

    p = sympy.Symbol('p')

    # The derived formula
    formula = (p**n - (1 - p)**n) / (2 * p - 1)

    # Use sympy to expand the formula into its polynomial form.
    # The division is exact, so we can use expand().
    poly = sympy.expand(formula)

    # The final equation is the polynomial expression for f_n(p).
    print(f"For n={n}, the final equation is:")
    # Use sympy.pretty_print for a nicely formatted output if in a compatible terminal
    # For this context, a string representation is better.
    poly_str = str(poly).replace('**', '^')
    print(f"f_{n}(p) = {poly_str}")

    # As requested, output each number (coefficient) in the final equation.
    poly_obj = sympy.Poly(poly, p)
    coeffs = poly_obj.all_coeffs()

    print("\nThe numbers (coefficients from highest power to lowest) in this equation are:")
    for c in coeffs:
        print(c)

if __name__ == '__main__':
    # The user did not specify a value for 'n'.
    # We will demonstrate the function for n=19, a value not in the provided list.
    generate_and_print_polynomial(19)
