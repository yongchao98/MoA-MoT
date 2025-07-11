import sympy

def poly_to_equation(p, var='X'):
    """Helper function to format a sympy polynomial into a string equation."""
    eq_parts = []
    # .all_coeffs() returns coefficients from highest degree to lowest
    coeffs = p.all_coeffs()
    degree = len(coeffs) - 1
    for i, coeff in enumerate(coeffs):
        power = degree - i
        if coeff == 0:
            continue
        # Format the term
        if power > 1:
            term = f"({coeff}) * {var}**{power}"
        elif power == 1:
            term = f"({coeff}) * {var}"
        else: # power == 0
            term = f"({coeff})"
        eq_parts.append(term)
    return " + ".join(eq_parts) + " = 0"

def main():
    # Let's define the symbolic variable for the polynomial
    X = sympy.Symbol('X')
    
    # Consider n=3. A matrix is derogatory if deg(pi_M) < 3.
    # The map M -> pi_M is discontinuous at derogatory matrices.
    # Let's create an example of such a matrix.
    M0 = sympy.Matrix([[5, 1, 0], 
                       [0, 5, 0], 
                       [0, 0, 5]])

    # Calculate the minimal polynomial of M0
    pi_M0 = M0.minpoly(X)
    deg_M0 = sympy.degree(pi_M0, X)

    print(f"Let's test the continuity at a derogatory matrix M0 (n=3).")
    print(f"M0 = \n{M0}\n")
    print(f"The minimal polynomial of M0 is pi_M0(X).")
    print(f"The equation pi_M0(X) = 0 is: {poly_to_equation(pi_M0)}")
    print(f"The degree of pi_M0 is {deg_M0}, which is less than n=3.")
    print("This means M0 is a derogatory matrix. The theory predicts the map is discontinuous here.\n")

    print("To show this, we construct a sequence of matrices Mk converging to M0.")
    # Mk is a small perturbation of M0
    # Mk = [[5, 1,    0],
    #       [0, 5, 1/k],
    #       [0, 0,    5]]
    # Let's pick a large k to make Mk very close to M0
    k_val = 10000
    Mk = sympy.Matrix([[5, 1, 0], 
                       [0, 5, sympy.Rational(1, k_val)], 
                       [0, 0, 5]])
    
    print(f"For k = {k_val}, the matrix Mk is:\n{Mk.evalf(6)}\n")
    print(f"As k -> infinity, Mk converges to M0.")

    # Calculate the minimal polynomial of Mk
    pi_Mk = Mk.minpoly(X)
    deg_Mk = sympy.degree(pi_Mk, X)
    
    print("Let's look at the minimal polynomial of Mk, pi_Mk(X).")
    print(f"The equation pi_Mk(X) = 0 is: {poly_to_equation(pi_Mk)}")
    print(f"The degree of pi_Mk is {deg_Mk}.\n")

    print("Summary of the result:")
    print("The sequence of polynomials pi_Mk is constant for any k > 0.")
    print("So, the limit of the sequence is lim (k->inf) pi_Mk(X):")
    print(f"   {poly_to_equation(pi_Mk)}")
    print("\nBut the minimal polynomial of the limit matrix M0 is pi_M0(X):")
    print(f"   {poly_to_equation(pi_M0)}")
    
    print("\nSince lim (k->inf) pi_Mk(X) is not equal to pi_M0(X), the map is discontinuous at M0.")

if __name__ == "__main__":
    main()
