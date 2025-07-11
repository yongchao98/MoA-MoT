import sympy
from sympy import I, pi, exp, gamma, sstr, Symbol

def solve_infinite_product():
    """
    This function calculates the symbolic expression for the infinite product
    Product_{n=3 to oo} (1 - z^3/n^3) and prints the resulting formula.
    """
    # Define the symbolic variable z
    z = Symbol('z')

    # The infinite product from n=1 to infinity is related to the Gamma function.
    # The formula is: Prod_{n=1 to inf} (1 - z^3/n^3) = 1 / P_gamma
    # where P_gamma = gamma(1-z) * gamma(1-omega*z) * gamma(1-omega^2*z)
    # and omega = exp(I*2*pi/3) is a complex cube root of unity.
    omega = exp(I * 2 * pi / 3)
    gamma_product = gamma(1 - z) * gamma(1 - omega * z) * gamma(1 - omega**2 * z)

    # The product in the question starts from n=3.
    # We must divide the full product by the terms for n=1 and n=2.
    n_1 = 1
    term_n1 = (1 - z**3 / n_1**3)

    n_2 = 2
    term_n2 = (1 - z**3 / n_2**3)

    # The final expression is 1 divided by the product of the gamma terms and the first two terms.
    final_expression = 1 / (term_n1 * term_n2 * gamma_product)

    # The prompt asks to output each number in the final equation.
    # We will print the final formula, showing the numbers 1 and 2 from the excluded terms.
    # The starting index 3 is from the problem statement.
    
    print("The final expression for the infinite product is:")
    # We use sstr for a compact string representation of the sympy expression.
    # The numbers 1 and 2 are explicitly shown in the denominator.
    print(f"1 / ( (1 - z**3/{n_1}**3) * (1 - z**3/{n_2}**3) * {sstr(gamma_product)} )")

# Execute the function to print the result
solve_infinite_product()