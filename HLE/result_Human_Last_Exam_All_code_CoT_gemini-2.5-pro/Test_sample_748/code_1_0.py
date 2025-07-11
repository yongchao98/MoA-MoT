import numpy as np
import sympy

def demonstrate_discontinuity():
    """
    This function demonstrates the discontinuity of the map M -> pi_M at a derogatory matrix.
    We use the 2x2 zero matrix M_0, which is derogatory.
    We construct a sequence of matrices M_k that converges to M_0, but for which
    the limit of their minimal polynomials does not equal the minimal polynomial of M_0.
    """
    # The variable for the polynomials
    x = sympy.Symbol('x')
    n = 2  # The dimension of the matrix

    # M_0: The 2x2 zero matrix. It is derogatory.
    # Its characteristic polynomial is x^2. Its minimal polynomial is x.
    M0_np = np.array([[0.0, 0.0], [0.0, 0.0]])
    M0 = sympy.Matrix(M0_np)
    pi_M0 = M0.minimal_polynomial(x)

    print("--- Demonstration of Discontinuity ---")
    print(f"Let's test continuity at the derogatory matrix M_0 =\n{M0_np}")
    print(f"The minimal polynomial of M_0 is pi_M0(x) = {pi_M0}")
    print("-" * 35)

    # M_k: A sequence of matrices converging to M_0.
    # M_k = [[1/k, 0], [0, -1/k]]. As k -> infinity, M_k -> M_0.
    # For any k > 0, M_k has distinct eigenvalues 1/k and -1/k.
    # Thus, its minimal polynomial is (x - 1/k)(x + 1/k) = x^2 - 1/k^2.
    k = 1e6
    Mk_np = np.array([[1/k, 0.0], [0.0, -1/k]])
    Mk = sympy.Matrix(Mk_np)
    pi_Mk = Mk.minimal_polynomial(x)
    
    # As k -> infinity, pi_Mk -> x^2.
    limit_poly = x**2

    print(f"Consider the sequence M_k = [[1/k, 0], [0, -1/k]], which converges to M_0.")
    print(f"For a large k, e.g., k = {k:.0e}, M_k is:\n{Mk_np}")
    print(f"The minimal polynomial of M_k is pi_Mk(x) = {sympy.expand(pi_Mk)}")
    print(f"As k -> infinity, the limit of these polynomials is lim(pi_Mk) = {limit_poly}")
    print("-" * 35)

    print("For continuity, we need lim(pi_Mk) to be equal to pi_M0.")
    print(f"We found lim(pi_Mk) = {limit_poly}")
    print(f"We know pi_M0 = {pi_M0}")
    print(f"Clearly, {limit_poly} != {pi_M0}.")
    print("\nConclusion: The map is discontinuous at M_0.")
    print("\nTo be explicit, we compare their coefficients in the polynomial space P_2(x).")

    # Get coefficients as vectors. The standard basis is {x^2, x^1, x^0}.
    # We pad with zeros for lower-degree polynomials to represent them in the space P_n[X].
    coeffs_limit = sympy.Poly(limit_poly, x).all_coeffs()
    padded_coeffs_limit = [float(c) for c in ([0]*(n + 1 - len(coeffs_limit)) + coeffs_limit)]

    coeffs_M0 = sympy.Poly(pi_M0, x).all_coeffs()
    padded_coeffs_M0 = [float(c) for c in ([0]*(n + 1 - len(coeffs_M0)) + coeffs_M0)]

    print(f"\nThe 'final equation' we are checking is:  {limit_poly} = {pi_M0}")
    print("In vector form, this is comparing the coefficients:")
    print(f"Coefficients of '{limit_poly}': {padded_coeffs_limit}")
    print(f"Coefficients of '{pi_M0}':    {padded_coeffs_M0}")
    print("Since the coefficient vectors are not equal, the polynomials are not equal.")

if __name__ == '__main__':
    demonstrate_discontinuity()
