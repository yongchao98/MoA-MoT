import numpy as np
import sympy as sp

def get_min_poly_obj(matrix):
    """Computes the minimal polynomial of a matrix as a SymPy Poly object."""
    X = sp.Symbol('X')
    sympy_matrix = sp.Matrix(matrix)
    min_poly_expr = sympy_matrix.minpoly(X)
    return sp.Poly(min_poly_expr, X)

def print_poly_equation(poly_obj):
    """Prints a polynomial in a readable equation format."""
    X = poly_obj.gen
    equation = sp.Eq(poly_obj.as_expr(), 0)
    sp.pprint(equation, use_unicode=True)

# Step 1: Define a derogatory matrix M0 for n=3.
# This matrix has eigenvalues {2, 2, 3}. The eigenvalue 2 has geometric multiplicity 2,
# which makes the matrix derogatory.
M0 = np.array([[2.0, 0.0, 0.0],
               [0.0, 2.0, 0.0],
               [0.0, 0.0, 3.0]])

# Step 2: Compute and display the minimal polynomial of M0.
pi_M0_obj = get_min_poly_obj(M0)

print("--- Demonstration of Discontinuity ---")
print("Consider the derogatory matrix M0:")
sp.pprint(sp.Matrix(M0), use_unicode=True)

print("\nIts minimal polynomial, pi_M0(X), is:")
# The minimal polynomial is (X-2)(X-3) = X^2 - 5X + 6.
print_poly_equation(pi_M0_obj)
print(f"The degree of pi_M0 is {pi_M0_obj.degree()}, which is less than n=3.")
coeffs_M0 = pi_M0_obj.all_coeffs()
print(f"The coefficients of the polynomial (from highest degree) are: {coeffs_M0}\n")

# Step 3: Define a sequence of matrices Mk that converges to M0.
# We create Mk by adding a small term 1/k that couples the two '2' eigenvalues.
# This perturbation makes Mk non-derogatory for any k > 0.
print("Now, consider a sequence of matrices Mk = [[2, 1/k, 0], [0, 2, 0], [0, 0, 3]] which converges to M0.")
print("For any k > 0, Mk is non-derogatory.")

# The minimal polynomial of Mk is the same for all k > 0.
# It equals the characteristic polynomial, (X-2)^2 * (X-3).
Mk_sample = np.array([[2.0, 0.01, 0.0],  # Using a small value for 1/k
                      [0.0, 2.0, 0.0],
                      [0.0, 0.0, 3.0]])
pi_Mk_obj = get_min_poly_obj(Mk_sample)
coeffs_Mk = pi_Mk_obj.all_coeffs()

print("\nThe minimal polynomial for any Mk (k>0), pi_Mk(X), is:")
print_poly_equation(pi_Mk_obj)
print(f"The degree of pi_Mk is {pi_Mk_obj.degree()}, which is equal to n=3.")
print(f"The coefficients of the polynomial are: {coeffs_Mk}\n")

# Step 4: Show that the limit of polynomials does not match.
print("As k -> infinity, Mk converges to M0.")
print("The limit of the minimal polynomials, lim(k->inf) pi_Mk(X), is the polynomial shown just above.")
print("This limit polynomial is NOT equal to pi_M0(X).")
print("\nFinal Polynomial Equation for pi_M0(X):")
print(f"1*X**2 + ({coeffs_M0[1]})*X**1 + ({coeffs_M0[2]})*X**0 = 0")
print("\nFinal Polynomial Equation for lim(k->inf) pi_Mk(X):")
print(f"1*X**3 + ({coeffs_Mk[1]})*X**2 + ({coeffs_Mk[2]})*X**1 + ({coeffs_Mk[3]})*X**0 = 0")
print("\nThis demonstrates that the map M -> pi_M is discontinuous at the derogatory matrix M0.")
