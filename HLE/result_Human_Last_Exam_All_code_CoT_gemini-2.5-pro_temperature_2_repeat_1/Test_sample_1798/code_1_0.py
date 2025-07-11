import numpy as np
from scipy.integrate import quad

# Step 1: Find the equation of the circle passing through A(1,0), B(10,0), C(4,3).
# The center (h,k) lies on the perpendicular bisector of AB.
# The x-coordinate of the center h is the midpoint of the x-coordinates of A and B.
h = (1 + 10) / 2

# Substitute h and the coordinates of A(1,0) and C(4,3) into the circle equation
# (x-h)^2 + (y-k)^2 = r^2 to solve for k and r^2.
# From A: (1 - 5.5)^2 + (0 - k)^2 = r^2  => 20.25 + k^2 = r^2
# From C: (4 - 5.5)^2 + (3 - k)^2 = r^2  => 2.25 + (3-k)^2 = r^2
# Equating the expressions for r^2:
# 20.25 + k^2 = 2.25 + 9 - 6k + k^2
# 20.25 = 11.25 - 6k => 9 = -6k => k = -1.5
k = -1.5

# Substitute k back to find r^2:
# r_sq = 20.25 + (-1.5)^2
r_sq = 20.25 + 2.25

# Step 2: Define the function y = f(x).
# From (y-k)^2 = r^2 - (x-h)^2, we get y = k +/- sqrt(r^2 - (x-h)^2).
# We choose the '+' to ensure the circle passes through C(4,3).
def circle_func(x):
    """Represents the upper part of the circle y=f(x)."""
    val_under_sqrt = r_sq - (x - h)**2
    # This check prevents math domain errors for inputs outside the circle's x-range.
    if val_under_sqrt < 0:
        return 0.0
    return k + np.sqrt(val_under_sqrt)

# Step 3: Find alpha. The total area under the PDF must be 1.
# So, alpha * integral of f(x) from 1 to 10 = 1.
total_integral, _ = quad(circle_func, 1, 10)
alpha = 1 / total_integral

# Step 4: Calculate P(X < 3).
# P(X < 3) is the integral of the PDF from 1 to 3.
integral_lt_3, _ = quad(circle_func, 1, 3)
prob_lt_3 = alpha * integral_lt_3

# As requested, output the numbers in the final equations.
print("Equation of the circle: (x - h)^2 + (y - k)^2 = r^2")
print(f"(x - {h})^2 + (y - ({k}))^2 = {r_sq}")

print("\nValue of the normalization constant alpha:")
print(f"alpha = {alpha}")

print("\nValue of the probability P(X < 3):")
print(f"P(X < 3) = {prob_lt_3}")