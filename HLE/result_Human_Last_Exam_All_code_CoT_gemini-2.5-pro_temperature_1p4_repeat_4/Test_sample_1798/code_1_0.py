import numpy as np
from scipy.integrate import quad

# Step 1: Define the parameters of the circle from the problem description.
# The points are A=(1,0), B=(10,0), and C=(4,3).
# The perpendicular bisector of the chord AB gives the x-coordinate of the center:
# h = (1 + 10) / 2 = 5.5.
# Using the general circle equation (x-h)^2 + (y-k)^2 = r^2 with points A and C:
# For A(1,0): (1 - 5.5)^2 + (0 - k)^2 = r^2 => 20.25 + k^2 = r^2
# For C(4,3): (4 - 5.5)^2 + (3 - k)^2 = r^2 => 2.25 + (3 - k)^2 = r^2
# Equating the expressions for r^2: 20.25 + k^2 = 2.25 + 9 - 6k + k^2
# This simplifies to 9 = -6k, so k = -1.5.
# The radius squared r^2 = 20.25 + (-1.5)^2 = 20.25 + 2.25 = 22.5.
h = 5.5
k = -1.5
r_squared = 22.5

# Step 2: Define the function y = f(x).
# We solve (y - k)^2 = r_squared - (x - h)^2 for y.
# y = k +/- sqrt(r_squared - (x - h)^2)
# y = -1.5 +/- sqrt(22.5 - (x - 5.5)^2)
# To pass through C(4,3), we must choose the '+' sign:
# -1.5 + sqrt(22.5 - (4-5.5)^2) = -1.5 + sqrt(20.25) = -1.5 + 4.5 = 3.
def f(x):
    """
    Represents the upper part of the circle equation passing through the given points.
    """
    discriminant = r_squared - (x - h)**2
    # The domain [1, 10] is valid for the square root.
    return k + np.sqrt(discriminant)

# Step 3: Calculate alpha.
# For d_X(x) = alpha * f(x) to be a PDF on [1, 10], its integral over this domain must be 1.
# alpha * integral(f(x) from 1 to 10) = 1 => alpha = 1 / integral(f(x) from 1 to 10).
total_integral, _ = quad(f, 1, 10)
alpha = 1 / total_integral

# Step 4: Calculate P(X < 3).
# This is the integral of the PDF from the start of the domain (1) to 3.
# P(X < 3) = integral(alpha * f(x) from 1 to 3).
partial_integral, _ = quad(f, 1, 3)
prob_X_lt_3 = alpha * partial_integral

# Step 5: Print the results.
print("The equation of the circle is (x - h)^2 + (y - k)^2 = r^2.")
print("The numbers in this equation are:")
print(f"h = {h}")
print(f"k = {k}")
print(f"r^2 = {r_squared}")
print(f"\nThe equation can be written as: (x - {h})^2 + (y - ({k}))^2 = {r_squared}")
print("\nThe values for the probability distribution are:")
print(f"alpha = {alpha}")
print(f"P(X < 3) = {prob_X_lt_3}")