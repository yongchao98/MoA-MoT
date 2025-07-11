import math
from scipy.integrate import quad

# Step 1: Determine the parameters of the circle passing through A(1,0), B(10,0), and C(4,3).
# The general equation is (x-h)^2 + (y-k)^2 = r^2.
# The x-coordinate of the center (h) is on the perpendicular bisector of AB.
h = (1 + 10) / 2

# Using points A(1,0) and C(4,3) and h=5.5, we solve for k and r^2.
# (1-5.5)^2 + (0-k)^2 = r^2  => 20.25 + k^2 = r^2
# (4-5.5)^2 + (3-k)^2 = r^2  => 2.25 + 9 - 6k + k^2 = r^2
# Equating them: 20.25 + k^2 = 11.25 - 6k + k^2 => 9 = -6k => k = -1.5
k = -1.5

# Substitute h and k back into the equation for point A to find r^2.
r_squared = (1 - h)**2 + (0 - k)**2

# Step 2: Define the function y = f(x).
# This is the upper arc of the circle, so we solve for y and take the positive root.
# (y-k)^2 = r^2 - (x-h)^2 => y = k + sqrt(r^2 - (x-h)^2)
def f(x):
    """Equation of the upper arc of the circle."""
    radicand = r_squared - (x - h)**2
    # The domain [1, 10] ensures the radicand is non-negative.
    return k + math.sqrt(radicand)

# Step 3: Calculate alpha.
# The PDF d_X(x) = alpha * f(x), and its integral over [1, 10] must be 1.
# alpha * integral(f(x) from 1 to 10) = 1 => alpha = 1 / integral(f(x) from 1 to 10)
integral_total, _ = quad(f, 1, 10)
alpha = 1 / integral_total

# Step 4: Calculate P(X < 3).
# P(X < 3) = integral(d_X(x) from 1 to 3) = alpha * integral(f(x) from 1 to 3).
integral_for_prob, _ = quad(f, 1, 3)
probability_X_less_than_3 = alpha * integral_for_prob

# Step 5: Print the results, showing the numbers in the final calculation.
print(f"The equation of the circle is (x - {h})^2 + (y - ({k}))^2 = {r_squared:.4f}")
print(f"The function y=f(x) is y = {k} + sqrt({r_squared:.4f} - (x - {h})^2)")
print("-" * 30)
print(f"To find alpha, we compute the integral of f(x) from 1 to 10:")
print(f"Integral_total = {integral_total:.4f}")
print(f"alpha = 1 / Integral_total = 1 / {integral_total:.4f} = {alpha:.4f}")
print("-" * 30)
print(f"To find P(X < 3), we compute the integral of f(x) from 1 to 3:")
print(f"Integral_prob = {integral_for_prob:.4f}")
print(f"P(X < 3) = alpha * Integral_prob = {alpha:.4f} * {integral_for_prob:.4f} = {probability_X_less_than_3:.4f}")

# Final answer in the required format
# <<<alpha, P(X<3)>>>
alpha_val = alpha
prob_val = probability_X_less_than_3