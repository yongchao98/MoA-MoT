import sympy

def illustrate_discontinuity():
    """
    This function demonstrates the discontinuity of the map M -> pi_M
    at a derogatory matrix.
    """
    X = sympy.Symbol('X')

    # Consider a derogatory matrix M, the 2x2 identity matrix.
    # It is derogatory because its minimal polynomial has degree 1 < 2.
    M = sympy.Matrix([[1, 0], [0, 1]])
    pi_M_poly = sympy.Poly(M.minimal_polynomial(X))
    chi_M_poly = sympy.Poly(M.charpoly(X))

    print("--- Analysis of a Derogatory Matrix M ---")
    print("Matrix M:")
    sympy.pprint(M)
    print(f"Minimal polynomial pi_M(X): {pi_M_poly.as_expr()}")
    print(f"Coefficients of pi_M(X) (descending powers): {pi_M_poly.all_coeffs()}")
    print(f"Characteristic polynomial chi_M(X): {chi_M_poly.as_expr()}")
    print(f"Coefficients of chi_M(X) (descending powers): {chi_M_poly.all_coeffs()}")
    print("-" * 45)

    # Consider a sequence of matrices M_eps converging to M as eps -> 0.
    # For any eps > 0, M_eps has distinct eigenvalues (1 and 1+eps),
    # so it is non-derogatory. Its minimal polynomial equals its characteristic polynomial.
    eps = sympy.Symbol('eps', real=True, positive=True)
    M_eps = sympy.Matrix([[1, 0], [0, 1 + eps]])
    pi_M_eps_poly = sympy.Poly(M_eps.minimal_polynomial(X))

    print("\n--- Analysis of a Nearby Non-Derogatory Matrix M_eps ---")
    print("Matrix M_eps (converges to M as eps -> 0):")
    sympy.pprint(M_eps)
    print(f"Minimal polynomial pi_M_eps(X): {pi_M_eps_poly.as_expr()}")
    print(f"Coefficients of pi_M_eps(X): {pi_M_eps_poly.all_coeffs()}")
    print("-" * 45)


    # Now, let's find the limit of the minimal polynomial of M_eps as eps -> 0.
    # We do this by taking the limit of each coefficient.
    limit_coeffs = [sympy.limit(c, eps, 0) for c in pi_M_eps_poly.all_coeffs()]
    limit_poly_expr = sympy.Poly.from_list(limit_coeffs, X).as_expr()

    print("\n--- The Limit ---")
    print(f"As eps -> 0, the coefficients of pi_M_eps(X) converge to: {limit_coeffs}")
    print(f"This corresponds to the polynomial: {limit_poly_expr}")
    print("\n--- Conclusion ---")
    print("The limit of the minimal polynomials of the sequence,")
    print(f"lim (eps->0) pi_M_eps(X) = {limit_poly_expr}")
    print("is equal to the characteristic polynomial of M, chi_M(X).")
    print("However, this is NOT equal to the minimal polynomial of M,")
    print(f"pi_M(X) = {pi_M_poly.as_expr()}.")
    print("\nThis demonstrates that the map M -> pi_M is discontinuous at M.")

if __name__ == '__main__':
    illustrate_discontinuity()