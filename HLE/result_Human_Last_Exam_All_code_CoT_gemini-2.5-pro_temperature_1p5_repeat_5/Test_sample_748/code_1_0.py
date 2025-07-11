import sympy as sp
import numpy as np

def print_poly_equation(poly, x, name):
    """
    Prints a polynomial equation with each coefficient clearly displayed.
    Example: pi_M(x) = 1.0*x^2 + -2.0*x + 1.0
    """
    poly_obj = sp.Poly(poly, x)
    coeffs = poly_obj.all_coeffs()
    degree = poly_obj.degree()
    
    terms = []
    for i, coeff in enumerate(coeffs):
        power = degree - i
        # Format the term string
        if power > 1:
            term = f"{float(coeff)}*x^{power}"
        elif power == 1:
            term = f"{float(coeff)}*x"
        else:
            term = f"{float(coeff)}"
        terms.append(term)
    
    equation_str = f"{name}(x) = " + " + ".join(terms)
    # A bit of formatting for cleaner output
    equation_str = equation_str.replace("+ -", "- ")
    print(equation_str)


# --- Main script ---
# Set the dimension of the matrices
n = 3

# Define a symbolic variable for the polynomial
x = sp.Symbol('x')

# 1. Define a derogatory matrix M0
# We choose the 3x3 identity matrix. It's derogatory because its minimal
# polynomial has degree 1, which is less than n=3.
M0 = sp.eye(n)

# Compute the minimal polynomial of M0
pi_M0 = M0.minpoly(x)

print("--- Analysis of a derogatory matrix M0 ---")
print(f"Let M0 be the {n}x{n} identity matrix:")
print(np.array(M0.tolist(), dtype=float))
print(f"\nThe minimal polynomial of M0 is:")
print_poly_equation(pi_M0, x, "pi_M0")

# 2. Construct a non-derogatory matrix Mk close to M0.
# We choose a small perturbation of M0. Let Mk = I + (1/k)*C where C is
# a nilpotent matrix (companion matrix for X^n).
# This matrix Mk will be non-derogatory.
C = sp.Matrix(np.diag(np.ones(n - 1), k=-1))

# Let's pick a large value for k to make the perturbation small
k = 1000
Mk = M0 + (sp.S(1)/k) * C

# Compute the minimal polynomial of Mk
pi_Mk = Mk.minpoly(x)

print("\n--- Analysis of a perturbed matrix Mk ---")
print(f"Let's construct a matrix Mk very close to M0 (with k={k}):")
print(np.array(Mk.tolist(), dtype=float))
print(f"\nThe minimal polynomial of this perturbed matrix Mk is:")
print_poly_equation(pi_Mk, x, "pi_Mk")


# 3. Conclusion
print("\n--- Conclusion ---")
print(f"The matrix Mk is very close to M0 (the difference is of order 1/k = {1.0/k}).")
print("However, their minimal polynomials are quite different:")
print(f"pi_M0(x) has degree {sp.degree(pi_M0)}.")
print(f"pi_Mk(x) has degree {sp.degree(pi_Mk)}.")
print("As k -> infinity, Mk -> M0, but lim pi_Mk(x) != pi_M0(x).")
print("This illustrates the discontinuity of the minimal polynomial map at the derogatory matrix M0.")
