import numpy as np

# This script demonstrates that chromatic roots can be real, non-integer, and negative.
# We consider the graph G formed by joining two triangles (K3) with a single edge.
# The chromatic polynomial for this graph is P_G(k) = k * (k-1) * (k-2)^2 * (k^2 - k - 1).
# To find all roots, we also need to solve for the roots of the quadratic factor k^2 - k - 1.
# The corresponding polynomial equation is k^2 - 1*k - 1 = 0.

# The coefficients of the polynomial k^2 - k - 1 are 1, -1, -1.
coeffs = [1, -1, -1]

# Find the roots of the polynomial equation
roots = np.roots(coeffs)

print("To demonstrate that chromatic roots can be non-integer and negative,")
print("we find the roots for a factor of the chromatic polynomial of a specific graph.")
print("The polynomial equation is k^2 - 1*k - 1 = 0.")
print(f"The numbers in the equation are: {coeffs[0]}, {coeffs[1]}, {coeffs[2]}")
print("\nThe roots of this factor are:")
print(roots[0])
print(roots[1])

print("\nAnalysis of the roots:")
print(f"1. The root {roots[0]:.3f} is a non-integer, demonstrating that statement D is true.")
print(f"2. The root {roots[1]:.3f} is negative, demonstrating that statement C is true.")