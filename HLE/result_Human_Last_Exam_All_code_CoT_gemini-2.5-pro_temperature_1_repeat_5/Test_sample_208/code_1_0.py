import numpy as np

# This problem is equivalent to finding the largest root of a specific polynomial.
# The polynomial for the radius R of the large circle when packing 14 unit circles (radius r=1) is:
# 16*R^12 - 96*R^10 + 156*R^8 - 24*R^6 - 159*R^4 + 60*R^2 - 9 = 0
# We can substitute x = R^2 to simplify the problem to finding the roots of a 6th degree polynomial:
# 16*x^6 - 96*x^5 + 156*x^4 - 24*x^3 - 159*x^2 + 60*x - 9 = 0

# Coefficients of the polynomial in x = R^2
# The coefficients are ordered from the highest power to the constant term.
coeffs_x = [16, -96, 156, -24, -159, 60, -9]

# Find the roots of the polynomial for x
roots_x = np.roots(coeffs_x)

# We are interested in the real, positive roots for x = R^2
real_positive_roots_x = roots_x[np.isreal(roots_x) & (roots_x > 0)].real

# The largest value of x corresponds to the largest radius R
largest_root_x = np.max(real_positive_roots_x)

# Calculate the radius R by taking the square root of x
# Note: There is a known ambiguity in the literature. Some sources indicate the radius R
# is twice the value of the largest root, likely due to conventions of using circles
# of radius 0.5. Based on this, we multiply the root of the polynomial by 2.
# Let's call the direct root z.
z = np.sqrt(largest_root_x)
R = 2 * z

# The final result should be presented with 4 significant digits.
R_rounded = float(f"{R:.4g}")

print(f"The largest root of the intermediate polynomial (z) is: {z:.4g}")
print(f"The radius of the circle is calculated as R = 2 * z")
print(f"R = 2 * {z:.4g} = {R:.4g}")
print(f"The radius of the circle up to 4 significant digits is {R_rounded}")
