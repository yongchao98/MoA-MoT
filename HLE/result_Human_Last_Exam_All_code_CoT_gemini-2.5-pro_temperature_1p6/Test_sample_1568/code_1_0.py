def solve_infinite_product():
    """
    This function provides the closed-form expression for the infinite product:
    Product_{n=3 to infinity} (1 - z^3/n^3)

    The solution is derived using the Weierstrass product representation for the
    Gamma function generalized for k-th roots of unity.
    
    1. The identity for the product starting from n=1 is:
       Product_{n=1 to inf} (1 - z^3/n^3) = 1 / (Gamma(1-z) * Gamma(1-z*w) * Gamma(1-z*w^2))
       where w = exp(i*2*pi/3).

    2. To get the product from n=3, we divide the expression above by the terms for n=1 and n=2.
       Term n=1: (1 - z^3/1^3)
       Term n=2: (1 - z^3/2^3)

    3. This script will print the resulting final expression.
    """
    
    # The final equation is of the form: Product = Expression.
    # The prompt provides the LHS, so we only need to provide the RHS.
    # The final expression below contains all the numbers of the final equation,
    # such as 1, 3, 8, 2, 4, etc., as requested.
    
    final_expression = "1 / ((1 - z**3) * (1 - z**3 / 8) * Gamma(1 - z) * Gamma(1 - z*exp(2*pi*I/3)) * Gamma(1 - z*exp(4*pi*I/3)))"
    
    print(final_expression)

solve_infinite_product()