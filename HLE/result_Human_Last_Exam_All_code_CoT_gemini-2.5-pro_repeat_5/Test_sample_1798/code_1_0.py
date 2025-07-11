import numpy as np
from scipy import integrate

# Step 1: Define the function f(x) based on the circle's equation.
# The equation of the circle is (x - h)^2 + (y - k)^2 = r^2
# with center (h,k) = (5.5, -1.5) and radius squared r_squared = 22.5.
# We solve for y, taking the positive root to pass through C(4,3).
# y = f(x) = k + sqrt(r_squared - (x - h)^2)
def f(x):
    """
    Represents the upper part of the circle passing through A, B, and C.
    """
    h = 5.5
    k = -1.5
    r_squared = 22.5
    # The term inside the square root is non-negative for x in [1, 10].
    sqrt_term = np.sqrt(r_squared - (x - h)**2)
    return k + sqrt_term

# Step 2: Find the normalization constant alpha.
# For d_X(x) = alpha * f(x) to be a PDF, its integral over [1, 10] must be 1.
# alpha * integral(f(x) from 1 to 10) = 1
# alpha = 1 / integral(f(x) from 1 to 10)

# Calculate the total integral of f(x) from 1 to 10 using numerical integration.
total_integral, _ = integrate.quad(f, 1, 10)

# Calculate alpha
alpha = 1 / total_integral

# Step 3: Calculate the probability P(X < 3).
# P(X < 3) = integral(alpha * f(x) from 1 to 3)
#           = alpha * integral(f(x) from 1 to 3)

# Calculate the partial integral of f(x) from 1 to 3.
partial_integral, _ = integrate.quad(f, 1, 3)

# Calculate the probability. This can also be seen as the ratio of the partial integral to the total integral.
probability_X_lt_3 = alpha * partial_integral

# Print the results and the numbers used in the final equations.
print("The equation for the upper semi-circle is y = f(x), where:")
print(f"f(x) = -1.5 + sqrt(22.5 - (x - 5.5)^2)")
print("\nTo find alpha, we solve the equation: alpha * Integral(f(x) dx from 1 to 10) = 1")
print(f"The value of the integral of f(x) from 1 to 10 is: {total_integral:.8f}")
print(f"So, alpha = 1 / {total_integral:.8f}")
print(f"alpha = {alpha:.8f}")

print("\nTo find P(X < 3), we solve the equation: P(X < 3) = alpha * Integral(f(x) dx from 1 to 3)")
print(f"The value of the integral of f(x) from 1 to 3 is: {partial_integral:.8f}")
print(f"So, P(X < 3) = {alpha:.8f} * {partial_integral:.8f}")
print(f"P(X < 3) = {probability_X_lt_3:.8f}")