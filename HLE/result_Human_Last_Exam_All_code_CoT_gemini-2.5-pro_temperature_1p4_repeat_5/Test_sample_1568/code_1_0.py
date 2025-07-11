import sympy
from sympy import I, pi, exp

def solve_infinite_product():
    """
    This function calculates the closed-form expression for the given infinite product
    and prints the result as a formatted equation.
    """
    # Define the symbolic variable z.
    z = sympy.Symbol('z')

    # Define the principal cube root of unity, omega.
    omega = exp(2 * pi * I / 3)

    # The left-hand side of the equation we want to solve.
    lhs_str = "Product_{n=3 to oo} (1 - z**3/n**3)"

    # The result is derived from the well-known product identity for the Gamma function.
    # The product from n=1 to infinity is 1 / (Gamma(1-z)*Gamma(1-z*w)*Gamma(1-z*w^2)).
    # We divide this by the terms for n=1 and n=2.
    # Term n=1: (1 - z^3/1^3) = 1 - z**3
    # Term n=2: (1 - z^3/2^3) = 1 - z**3/8 = (8 - z**3)/8
    # The denominator becomes (1 - z**3)*(8 - z**3)/8 * Product_of_Gammas
    # This simplifies to 8 / ((1-z**3)*(8-z**3)*Product_of_Gammas).

    # To satisfy the prompt's instruction to "output each number in the final equation",
    # we define the numerical components of the final expression here.
    numerator_val = 8
    
    # Components of the (1-z^3) term
    c1 = 1
    p1 = 3
    
    # Components of the (8-z^3) term
    c2 = 8
    p2 = 3

    # Construct the symbolic expression for the final result.
    gamma_product = sympy.gamma(1 - z) * sympy.gamma(1 - z * omega) * sympy.gamma(1 - z * omega**2)
    denominator_z_terms = (c1 - z**p1) * (c2 - z**p2)
    
    # The full expression. We use sympy.Mul and sympy.Pow for proper representation.
    # result = 8 / ( (1-z**3)*(8-z**3) * gamma(1-z)*gamma(1-z*omega)*gamma(1-z*omega**2) )
    result_expr = sympy.Mul(numerator_val, sympy.Pow(denominator_z_terms * gamma_product, -1))

    # Print the final equation in a readable format using sympy.pprint.
    print(f"{lhs_str} =")
    sympy.pprint(result_expr, use_unicode=True)

solve_infinite_product()