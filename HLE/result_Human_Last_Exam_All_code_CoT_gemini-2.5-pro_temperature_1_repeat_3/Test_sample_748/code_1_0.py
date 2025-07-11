import numpy as np
import sympy as sp

def get_minimal_polynomial(matrix, var='x'):
    """Computes the minimal polynomial of a matrix using sympy."""
    # sympy's minpoly function is reliable for this.
    M_sp = sp.Matrix(matrix)
    x = sp.Symbol(var)
    min_poly = M_sp.minpoly(x)
    return min_poly

def get_coeffs_vector(poly, max_deg, var='x'):
    """Extracts coefficients of a sympy polynomial into a list."""
    # We represent a polynomial as a vector of its coefficients.
    x = sp.Symbol(var)
    # The Poly class helps in extracting coefficients.
    p = sp.Poly(poly, x)
    # Initialize a zero vector for coefficients up to max_deg
    coeffs = [0] * (max_deg + 1)
    # .all_coeffs() gives pairs of (degree, coefficient)
    # We map them to the correct position in our vector.
    # Note: sympy returns highest degree first, so we reverse to get [c0, c1, ...]
    all_c = p.all_coeffs()
    current_deg = p.degree()
    for c in all_c:
        coeffs[current_deg] = c
        current_deg -=1
    return coeffs

# Set the dimension of the matrix space
n = 3
x_sym = sp.Symbol('x')

# --- Step 1: Define a derogatory matrix M0 ---
# A matrix is derogatory if its minimal polynomial has degree < n.
# The 3x3 zero matrix is a simple example.
M0 = np.zeros((n, n))
pi_M0 = get_minimal_polynomial(M0)
deg_pi_M0 = sp.degree(pi_M0, gen=x_sym)

print(f"Let n = {n}.")
print("We analyze the continuity at a derogatory matrix M0.")
print(f"Let M0 be the {n}x{n} zero matrix:\n{M0}")
print(f"The minimal polynomial of M0 is pi_M0(x) = {pi_M0}.")
print(f"The degree of pi_M0 is {deg_pi_M0}, which is less than n={n}. So M0 is derogatory.")
print("-" * 40)

# --- Step 2: Define a sequence of non-derogatory matrices Mk -> M0 ---
# A matrix is non-derogatory if its minimal polynomial has degree n.
# We choose a sequence of matrices Mk that are non-derogatory for any k,
# but converge to the derogatory matrix M0 as k -> infinity.
# A scaled nilpotent Jordan block is a good example.
k_val = 1000.0
Mk = np.array([
    [0, 1/k_val, 0],
    [0, 0,   1/k_val],
    [0, 0,   0]
], dtype=float)
# For any k > 0, the minimal polynomial of Mk is x^3.
pi_Mk = get_minimal_polynomial(Mk)
deg_pi_Mk = sp.degree(pi_Mk, gen=x_sym)

print("Consider a sequence of matrices Mk converging to M0 as k -> infinity.")
print(f"Let's take an element of this sequence, for k = {k_val}:")
print(f"Mk =\n{Mk}")
print(f"The minimal polynomial of this Mk is pi_Mk(x) = {pi_Mk}.")
print(f"The degree of pi_Mk is {deg_pi_Mk}, which is equal to n={n}. So Mk is non-derogatory.")
print("-" * 40)

# --- Step 3: Check for continuity ---
# The map is continuous if lim(pi_Mk) = pi(lim(Mk)).
# lim(Mk) = M0, so we check if lim(pi_Mk) = pi_M0.

# The minimal polynomial of Mk is x^3 for any k > 0.
# The sequence of polynomials pi_Mk is constant: x^3, x^3, x^3, ...
# The limit of this sequence is trivially x^3.
limit_pi_Mk = sp.sympify("x**3")

# We represent polynomials in C_n[x] by their coefficient vectors.
# For n=3, the space is C_3[x] with basis (1, x, x^2, x^3).
coeffs_pi_M0 = get_coeffs_vector(pi_M0, max_deg=n)
coeffs_limit_pi_Mk = get_coeffs_vector(limit_pi_Mk, max_deg=n)

print("To check for continuity, we compare the limit of the polynomials")
print("with the polynomial of the limit matrix.")
print("\nRepresenting polynomials by their coefficient vectors (for basis 1, x, x^2, x^3):")
# The final equation comparison:
print(f"pi( lim Mk ) = pi(M0) = {pi_M0}")
print(f"Coefficient vector: {coeffs_pi_M0}")

print(f"lim( pi_Mk ) = {limit_pi_Mk}")
print(f"Coefficient vector: {coeffs_limit_pi_Mk}")

print("\nSince the polynomials (and their coefficient vectors) are not equal,")
print("the map is discontinuous at the derogatory matrix M0.")
