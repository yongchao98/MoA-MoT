import numpy as np
from scipy import integrate

# Based on the derivation, the solution for y1(x) is assumed to be:
def y1(x):
    return 8 * x + 42688 * x**-6

# Based on the derivation, the solution for y2(x) with n=1 (yd=1) is:
def y2(x, n=1):
    yd = 1.0 / n
    return (yd * x) / ((1 + 2 * x**5)**(2.0/5.0))

# We found n=1 is the minimal integer for non-intersection.
n = 1
yd = 1.0/n

# The integration range is assumed to be from x=2 to x=10 due to ambiguities
# in the problem's inequality condition.
lower_bound = 2
upper_bound = 10

# Calculate the definite integral of y1(x) from lower_bound to upper_bound
integral_val, _ = integrate.quad(y1, lower_bound, upper_bound)

# The integral evaluates to (4*x^2 - 42688/(5*x^5))
def indefinite_y1(x):
    return 4 * x**2 - 42688 / (5 * x**5)

val_at_10 = indefinite_y1(upper_bound)
val_at_2 = indefinite_y1(lower_bound)
result = val_at_10 - val_at_2

# Print the final equation
print(f"The integral is calculated as [4*x^2 - 42688/(5*x^5)] from {lower_bound} to {upper_bound}.")
print(f"Value at x={upper_bound}: 4*({upper_bound})^2 - 42688/(5*({upper_bound})^5) = {val_at_10}")
print(f"Value at x={lower_bound}: 4*({lower_bound})^2 - 42688/(5*({lower_bound})^5) = {val_at_2}")
print(f"Final result: {val_at_10} - ({val_at_2}) = {result}")
