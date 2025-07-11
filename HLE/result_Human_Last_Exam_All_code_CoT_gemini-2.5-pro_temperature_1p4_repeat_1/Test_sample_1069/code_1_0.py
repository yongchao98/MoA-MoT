import math

# The specific solution to the boundary-value problem is x(t) = (2 - exp(-t)) / cosh(t).
# We need to find the value of x(t) for t = ln(10^20).

# Define the terms based on t = ln(10^20)
# exp(t) = exp(ln(10^20)) = 10^20
e_t = 10.0**20

# exp(-t) = 1 / exp(t)
e_neg_t = 1.0 / e_t

# cosh(t) = (e^t + e^-t) / 2
cosh_t = (e_t + e_neg_t) / 2.0

# Numerator of the solution is 2 - e^-t
numerator = 2.0 - e_neg_t

# Calculate x(t)
x_t = numerator / cosh_t

print("The specific solution is x(t) = (2 - exp(-t)) / cosh(t).")
print("For t = ln(10^20), we calculate the following values:")

print(f"\nexp(-t) = {e_neg_t}")
print(f"cosh(t) = {cosh_t}")
print(f"Numerator (2 - exp(-t)) = {numerator}")

print(f"\nThe final calculation is:")
print(f"x(ln(10^20)) = {numerator} / {cosh_t}")

# Due to floating-point precision, the result is an approximation,
# but it is extremely close to the true value.
print(f"\nThe final result is:")
print(x_t)