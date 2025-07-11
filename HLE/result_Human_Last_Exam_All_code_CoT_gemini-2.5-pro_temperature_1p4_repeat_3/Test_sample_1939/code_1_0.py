import numpy as np

# This script demonstrates the solution for the case d=3.
# We first determine the Ehrhart polynomial p(z) for the polytope when d=3.
# p(z) is a polynomial of degree 3: p(z) = c3*z^3 + c2*z^2 + c1*z + c0.

# From the problem analysis, we have the following properties for d=3:
# 1. p(0) = 1, which implies c0 = 1.
# 2. p(1) = 2*d = 6.
# 3. p(-1) = 0.
# 4. The leading coefficient, c3, is the volume of the polytope, which is 2/3 for d=3.
c3 = 2/3
c0 = 1

# We have a system of linear equations for c1 and c2 from p(1) and p(-1):
# From p(1)=6: c3*1^3 + c2*1^2 + c1*1 + c0 = 6  => c3 + c2 + c1 + c0 = 6
# From p(-1)=0: c3*(-1)^3 + c2*(-1)^2 + c1*(-1) + c0 = 0 => -c3 + c2 - c1 + c0 = 0

# Substitute known values:
# 2/3 + c2 + c1 + 1 = 6   => c2 + c1 = 13/3
# -2/3 + c2 - c1 + 1 = 0  => c2 - c1 = -1/3

# We can solve this system.
# Adding the two equations: 2*c2 = 12/3 = 4 => c2 = 2.
# Substituting c2 back: 2 + c1 = 13/3 => c1 = 7/3.
c2 = 2
c1 = 7/3

# So the polynomial is p(z) = (2/3)z^3 + 2z^2 + (7/3)z + 1.
# The roots of p(z) are the same as the roots of 3*p(z) = 2z^3 + 6z^2 + 7z + 3.
coeffs = [c3, c2, c1, c0]
scaled_coeffs = [2, 6, 7, 3]

# We use numpy to find the roots of the polynomial.
roots = np.roots(scaled_coeffs)

print(f"The dimension is d=3.")
print(f"The coefficients of the Ehrhart polynomial p(z) are:")
print(f"c3 (leading coefficient) = {coeffs[0]:.4f}")
print(f"c2 = {coeffs[1]:.4f}")
print(f"c1 = {coeffs[2]:.4f}")
print(f"c0 = {coeffs[3]:.4f}")
print("\nSo, p(z) = (2/3)z^3 + 2z^2 + (7/3)z + 1")
print("The roots of this polynomial are:")
for r in roots:
    print(f"{r.real:.4f} + {r.imag:.4f}i")

# Check the real part of each root.
all_real_part_neg_one = all(np.isclose(r.real, -1) for r in roots)

print(f"\nDo all roots have a real part of -1? {all_real_part_neg_one}")
print("This confirms that for d=3, every root of p(z) has a real part of -1.")
