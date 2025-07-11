import numpy as np
from scipy.integrate import quad

# Step 1: Find the equation of the circle
# Given points A=(1,0), B=(10,0), C=(4,3)

# The perpendicular bisector of the line segment AB is a vertical line at x = (1+10)/2.
# So, the x-coordinate of the circle's center (h,k) is h = 5.5.
h = 5.5

# We use the fact that the distance from the center to any point on the circle is the radius r.
# r^2 = (x - h)^2 + (y - k)^2
# Using point A(1,0): r^2 = (1 - 5.5)^2 + (0 - k)^2 = 20.25 + k^2
# Using point C(4,3): r^2 = (4 - 5.5)^2 + (3 - k)^2 = 2.25 + (3 - k)^2
# Equating the two expressions for r^2:
# 20.25 + k^2 = 2.25 + 9 - 6k + k^2
# 20.25 = 11.25 - 6k
# 9 = -6k  => k = -1.5
k = -1.5

# Now we find r^2 using the value of k:
# r^2 = 20.25 + (-1.5)^2 = 20.25 + 2.25 = 22.5
r_squared = 22.5

print("Step 1: Finding the equation of the circle")
print(f"The center of the circle is (h, k) = ({h}, {k})")
print(f"The radius squared is r^2 = {r_squared}")
print(f"The equation of the circle is (x - {h})^2 + (y - ({k}))^2 = {r_squared}")
print(f"This is equivalent to: (x - 5.5)^2 + (y + 1.5)^2 = 22.5\n")

# Step 2: Define the function y = f(x)
# From the circle's equation: (y + 1.5)^2 = 22.5 - (x - 5.5)^2
# y = -1.5 +/- sqrt(22.5 - (x - 5.5)^2)
# We must use the '+' sign for the circle to pass through C=(4,3).
def f(x):
    """Represents the upper part of the circle y = f(x)."""
    return k + np.sqrt(r_squared - (x - h)**2)

print("Step 2: Defining the function y = f(x)")
print(f"The function is f(x) = {k} + sqrt({r_squared} - (x - {h})^2)\n")

# Step 3: Calculate alpha
# The probability density function d_X(x) = alpha * f(x) must integrate to 1 over its domain [1, 10].
# alpha * integral(f(x) from 1 to 10) = 1
# alpha = 1 / integral(f(x) from 1 to 10)

# We use numerical integration for accuracy.
total_integral, integral_error_total = quad(f, 1, 10)
alpha = 1 / total_integral

print("Step 3: Calculating alpha")
print(f"The definite integral of f(x) from 1 to 10 is: {total_integral}")
print(f"alpha = 1 / (integral of f(x))")
print(f"The value of alpha is: {alpha}\n")

# Step 4: Calculate the probability P(X < 3)
# P(X < 3) is the integral of the PDF from the start of the domain (1) to 3.
# P(X < 3) = integral(alpha * f(x) from 1 to 3) = alpha * integral(f(x) from 1 to 3)
partial_integral, integral_error_partial = quad(f, 1, 3)
probability = alpha * partial_integral

print("Step 4: Calculating P(X < 3)")
print(f"The definite integral of f(x) from 1 to 3 is: {partial_integral}")
print(f"P(X < 3) = alpha * (integral of f(x) from 1 to 3)")
print(f"The value of P(X < 3) is: {probability}\n")

print("--- Final Answer ---")
print(f"The value of alpha is: {alpha}")
print(f"The value of P(X < 3) is: {probability}")