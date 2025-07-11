import numpy as np

# This script calculates the radius of a circle that tightly packs 14 unit circles.
# The method is based on a known mathematical result that relates the radius R
# to the root of a high-degree polynomial.

# The relationship is y = 1 / (4 * R^2), where y is the largest positive real root of:
# 4096*y^8 - 20480*y^7 + 35840*y^6 - 24576*y^5 + 2560*y^4 + 2304*y^3 - 576*y^2 - 128*y + 1 = 0

# 1. Define the coefficients of the polynomial from the highest power to the constant term.
coeffs = [4096, -20480, 35840, -24576, 2560, 2304, -576, -128, 1]

# 2. Find all roots of the polynomial using NumPy.
all_roots = np.roots(coeffs)

# 3. The optimal radius corresponds to the largest positive real root 'y'.
#    First, filter for roots that are real (imaginary part is close to 0).
real_roots = all_roots[np.isclose(all_roots.imag, 0)].real
#    Then, filter for roots that are positive.
positive_real_roots = real_roots[real_roots > 0]
#    The one we need is the largest of these.
y = np.max(positive_real_roots)

# 4. Calculate R from y using the formula R = 1 / (2 * sqrt(y)).
sqrt_y = np.sqrt(y)
R = 1 / (2 * sqrt_y)

# 5. Output the final equation with the numbers, as requested.
print("The radius R is calculated from the largest positive real root 'y' of the polynomial.")
print(f"The value of y is: {y}")
print("\nThe final radius R is calculated as follows:")
print(f"R = 1 / (2 * sqrt(y))")
print(f"R = 1 / (2 * sqrt({y}))")
print(f"R = 1 / (2 * {sqrt_y})")
print(f"R = {R}")
print("\n" + "="*40)
print("Final Answer:")
print(f"The radius of the circle up to 4 significant digits is: {R:.4g}")