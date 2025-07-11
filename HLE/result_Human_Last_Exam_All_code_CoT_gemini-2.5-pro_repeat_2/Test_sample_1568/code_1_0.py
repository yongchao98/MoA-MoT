import sympy

def solve_infinite_product():
    """
    This function calculates the symbolic expression for the given infinite product
    using the Gamma function identity and prints the result.
    """
    # Define the symbolic variable 'z'.
    z = sympy.Symbol('z')

    # The problem is to evaluate the infinite product starting from n=3.
    # We use the known identity for the product from n=1.
    # The numbers in the equation are 1, 3, 2, 8, and the constants in Gamma and roots of unity.

    # Define the complex cubic roots of unity, w = exp(2*pi*I/3) and w^2 = exp(4*pi*I/3).
    # These contain the numbers 2, 3, 4, pi, and I.
    omega = sympy.exp(2 * sympy.pi * sympy.I / 3)
    omega_sq = sympy.exp(4 * sympy.pi * sympy.I / 3)

    # The product of Gamma functions in the denominator of the identity. The number 1 is in the argument.
    gamma_product = sympy.gamma(1 - z) * sympy.gamma(1 - z * omega) * sympy.gamma(1 - z * omega_sq)

    # The full infinite product starting from n=1.
    infinite_product_from_1 = 1 / gamma_product

    # The required product starts from n=3, so we divide by the terms for n=1 and n=2.
    # Term for n=1 contains the numbers 1 and 3: (1 - z**3 / 1**3).
    term_for_n_1 = (1 - z**3)
    # Term for n=2 contains the numbers 2 and 3, which result in 8: (1 - z**3 / 2**3).
    term_for_n_2 = (1 - z**3 / 8)

    # The final expression is the product from n=1 divided by the first two terms.
    final_expression = infinite_product_from_1 / (term_for_n_1 * term_for_n_2)

    # Print the final symbolic expression. The string representation will contain all the numbers.
    print(final_expression)

if __name__ == '__main__':
    solve_infinite_product()