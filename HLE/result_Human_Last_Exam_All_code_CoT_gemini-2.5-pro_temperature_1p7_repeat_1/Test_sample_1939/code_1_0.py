import numpy as np

# We found the Ehrhart polynomial for d=3 to be p(z) = (2z^3 + 6z^2 + 7z + 3)/3.
# The roots of this polynomial are the roots of the numerator.
# Let's find the roots of the cubic polynomial 2z^3 + 6z^2 + 7z + 3.

# Coefficients of the polynomial
coeffs = [2, 6, 7, 3]

# Find the roots
roots = np.roots(coeffs)

print("The Ehrhart polynomial for d=3 is p(z) = (2/3)z^3 + 2z^2 + (7/3)z + 1.")
print("The roots of the polynomial are:")
for root in roots:
    print(root)

print("\nThe real parts of the roots are:")
for root in roots:
    print(np.real(root))

print("\nAs can be seen, for d=3, every root has a real part of -1.")
print("For d=2, the polynomial is z^2 + 2z + 1 = (z+1)^2, whose only root is z=-1.")
print("This supports the conclusion that for any d, every root of the Ehrhart polynomial has a real part of -1.")
