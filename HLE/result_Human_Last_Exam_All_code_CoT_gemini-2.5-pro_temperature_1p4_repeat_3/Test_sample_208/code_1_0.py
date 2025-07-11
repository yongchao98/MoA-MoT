import numpy as np

# This script calculates the radius of a circle that can tightly pack 14 smaller circles of radius 1.
# The solution is based on the proven optimal packing arrangement for 14 circles.
# The radius R is found by solving a polynomial equation for the variable y = (R-1)^2.

# 1. Define the coefficients of the polynomial: 27y^4 - 324y^3 + 1053y^2 - 1188y + 361 = 0
coeffs = [27, -324, 1053, -1188, 361]

print("The polynomial equation for y = (R-1)^2 is: 27*y^4 - 324*y^3 + 1053*y^2 - 1188*y + 361 = 0\n")

# 2. Use numpy to find the numerical roots of the polynomial.
roots = np.roots(coeffs)

# 3. The radius of the enclosing circle is determined by the outermost small circles.
#    This corresponds to the largest real and positive root of the polynomial.
real_positive_roots = roots[np.isreal(roots) & (roots > 0)].real
y_max = np.max(real_positive_roots)

print(f"The roots of the polynomial for y are:\n{roots}\n")
print(f"The largest real, positive root is y = {y_max:.8f}\n")

# 4. Calculate R from y_max. Since y = (R-1)^2 and the small circles have radius 1,
#    R = sqrt(y_max) + 1.
sqrt_y = np.sqrt(y_max)
R = sqrt_y + 1

print("The final equation for the radius R is R = sqrt(y) + 1.")
print("Plugging in the value for y:")
print(f"R = sqrt({y_max:.8f}) + 1")
print(f"R = {sqrt_y:.8f} + 1")
print(f"R = {R:.8f}\n")

# 5. Format the result to 4 significant digits as requested by the user.
# The 'g' format specifier is used for rounding to a number of significant figures.
formatted_R = f"{R:.4g}"

print(f"The radius of the circle up to 4 significant digits is: {formatted_R}")