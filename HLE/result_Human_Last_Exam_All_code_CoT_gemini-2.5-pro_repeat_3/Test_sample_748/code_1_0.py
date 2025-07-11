import numpy as np
import sympy as sp

# Use sympy for polynomial manipulation
X = sp.Symbol('X')

def get_min_poly(matrix):
    """Computes the minimal polynomial of a matrix using sympy."""
    # Convert numpy matrix to sympy matrix
    M_sympy = sp.Matrix(matrix)
    # The minpoly function requires the variable for the polynomial
    min_poly = M_sympy.minpoly(X)
    return min_poly

# --- Step 1: Define a derogatory matrix M ---
# We choose the 2x2 identity matrix. It's derogatory because its
# minimal polynomial (X-1) has degree 1, which is less than n=2.
M = np.array([[1, 0], [0, 1]])

# --- Step 2: Calculate the minimal polynomial of M ---
pi_M = get_min_poly(M)

print("--- Discontinuity Demonstration ---")
print(f"Let M be the derogatory matrix:\n{M}\n")
print(f"The minimal polynomial of M is pi_M(X) = {sp.expand(pi_M)}\n")


# --- Step 3: Define a sequence of matrices M_k converging to M ---
# We choose M_k = [[1, 1/k], [0, 1]].
# As k -> infinity, M_k -> M.
# For any finite k > 0, M_k is not a scalar multiple of the identity,
# so it is non-derogatory. Its minimal polynomial is its characteristic polynomial (X-1)^2.
k = 10000 # Choose a large value for k to approximate the limit
M_k = np.array([[1, 1/k], [0, 1]])

# --- Step 4: Calculate the minimal polynomial of M_k ---
pi_M_k = get_min_poly(M_k)

print(f"Consider the sequence of non-derogatory matrices M_k = [[1, 1/k], [0, 1]].")
print(f"As k -> infinity, M_k converges to M.")
print(f"For a large k, let's take k = {k}, the matrix M_k is:\n{M_k}\n")
print(f"The minimal polynomial of M_k is pi_M_k(X) = {sp.expand(pi_M_k)}")

# --- Step 5: The limit of the minimal polynomials pi_M_k is (X-1)^2 ---
# For any k, pi_M_k = (X-1)^2. So the limit is clearly (X-1)^2.
limit_poly_expr = (X-1)**2
limit_poly = sp.expand(limit_poly_expr)

print(f"\nThe limit of pi_M_k(X) as k -> infinity is (X-1)^2 = {limit_poly}\n")


# --- Step 6: Compare the limit polynomial with pi_M ---
print("--- Conclusion ---")
print(f"The limit of the minimal polynomials of the sequence is: {limit_poly}")
print(f"The minimal polynomial of the limit matrix M is: {sp.expand(pi_M)}")
print(f"Since these are not equal, the map M -> pi_M is discontinuous at M.")

# --- Final Equation Comparison with Coefficients ---
# The prompt asks to output each number in the final equation.
# We interpret this as showing the coefficients of the polynomials.
limit_poly_coeffs = sp.Poly(limit_poly, X).all_coeffs()
pi_M_coeffs = sp.Poly(pi_M, X).all_coeffs()

print("\nFinal Equation Comparison:")
print(f"Limit polynomial equation: {limit_poly} = 0")
print("Coefficients (descending powers of X):", [float(c) for c in limit_poly_coeffs])
print(f"Minimal polynomial of M equation: {pi_M} = 0")
print("Coefficients (descending powers of X):", [float(c) for c in pi_M_coeffs])
