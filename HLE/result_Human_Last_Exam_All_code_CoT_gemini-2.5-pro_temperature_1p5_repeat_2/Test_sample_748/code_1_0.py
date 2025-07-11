import sympy

def demonstrate_discontinuity():
    """
    Demonstrates the discontinuity of the minimal polynomial map.
    We choose a 3x3 derogatory matrix M0 and a sequence of non-derogatory
    matrices Mk converging to M0. We then show that the minimal polynomial of Mk
    does not converge to the minimal polynomial of M0.
    """
    # Use symbolic variables for generality
    l, m, k, x = sympy.symbols('lambda mu k x')

    # Define a derogatory matrix M0.
    # It has eigenvalues lambda (with multiplicity 2) and mu.
    # The geometric multiplicity of lambda is 2, making it derogatory.
    M0 = sympy.Matrix([
        [l, 0, 0],
        [0, l, 0],
        [0, 0, m]
    ])
    
    # M0 is derogatory, so its minimal polynomial has degree less than 3.
    # We assume lambda != mu.
    pi_M0 = (x - l) * (x - m)
    pi_M0_expanded = sympy.expand(pi_M0)

    print("Chosen derogatory matrix M0:")
    sympy.pprint(M0)
    print("\nThe minimal polynomial of M0 is (x - lambda)(x - mu):")
    sympy.pprint(pi_M0)
    print("Expanded form: ")
    sympy.pprint(pi_M0_expanded)
    print("-" * 30)

    # Define a sequence of matrices Mk that converges to M0 as k -> infinity.
    # We add a small perturbation 1/k.
    Mk = sympy.Matrix([
        [l, 1/k, 0],
        [0, l,   0],
        [0, 0,   m]
    ])

    print("Sequence of matrices Mk converging to M0 (for large k):")
    sympy.pprint(Mk)

    # For any finite k, Mk is non-derogatory.
    # Its minimal polynomial is its characteristic polynomial.
    # The eigenvalues are still lambda, lambda, mu.
    char_poly_Mk = (x - l)**2 * (x - m)
    pi_Mk_expanded = sympy.expand(char_poly_Mk)
    
    print(f"\nThe minimal polynomial of Mk for any k > 0 is (x - lambda)^2 * (x - mu):")
    sympy.pprint(char_poly_Mk)
    print("Expanded form: ")
    sympy.pprint(pi_Mk_expanded)
    print("-" * 30)
    
    # As k approaches infinity, Mk approaches M0.
    # The minimal polynomial of Mk, however, approaches (x-l)^2(x-m).
    # This is not equal to the minimal polynomial of M0, which is (x-l)(x-m).
    print("Conclusion:")
    print("As k -> oo, Mk -> M0.")
    print("But lim (pi_Mk) != pi_M0:")
    print("Limit of coefficients of pi_Mk:", [1, -(2*l+m), l**2+2*l*m, -(l**2*m)])
    print("Coefficients of pi_M0:", [0, 1, -(l+m), l*m]) # Padded with 0 for degree 3
    print("\nThe polynomials are different, demonstrating the discontinuity at M0.")


if __name__ == '__main__':
    demonstrate_discontinuity()
