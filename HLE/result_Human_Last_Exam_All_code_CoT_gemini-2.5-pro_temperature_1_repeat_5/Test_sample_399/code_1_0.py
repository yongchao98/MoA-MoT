import math
from fractions import Fraction

# The initial value problem is x'(t) = (t-1) * x^2(t), with x(0) = -8.
# This is a separable differential equation. We can rewrite it as:
# dx / x^2 = (t - 1) dt

# Integrating both sides gives:
# -1/x = t^2/2 - t + C

# We use the initial condition x(0) = -8 to find the constant C.
t0 = 0
x0 = -8

# Substitute t0 and x0 into the integrated equation:
# -1/(-8) = 0^2/2 - 0 + C
# This simplifies to C = 1/8.
C = -1/x0 - (t0**2 / 2 - t0)
print(f"Step 1: Find the constant of integration, C.")
print(f"Using the initial condition x({t0}) = {x0}, we get C = -1/({x0}) = {C}")
print("-" * 50)

# The particular solution for x(t) is:
# -1/x = t^2/2 - t + 1/8
# Solving for x(t), we get:
# x(t) = 1 / (-t^2/2 + t - 1/8)

# Now, we need to find the value of x(1).
t_final = 1
print(f"Step 2: Use the particular solution to find x({t_final}).")
print(f"x(t) = 1 / (-t^2/2 + t - C)")
print(f"x({t_final}) = 1 / (-({t_final})^2/2 + {t_final} - {C})")
print("-" * 50)


# Define the numbers in the final equation to calculate x(1)
numerator = 1
term1_denom = -(t_final**2) / 2
term2_denom = t_final
term3_denom = -C

print("Step 3: Calculate the final value by plugging in the numbers.")
print("The final equation for x(1) is built from the following parts:")
print(f"Numerator: {numerator}")
print(f"Term 1 in denominator (-t^2/2): {term1_denom}")
print(f"Term 2 in denominator (+t): {term2_denom}")
print(f"Term 3 in denominator (-C): {term3_denom}")
print("-" * 50)

# Calculate the final result
denominator = term1_denom + term2_denom + term3_denom
x_final = numerator / denominator

print(f"The value of the denominator is: {term1_denom} + {term2_denom} + {term3_denom} = {denominator}")
print(f"The final value is x(1) = {numerator} / {denominator} = {x_final}")

# Display the exact fraction for clarity
final_fraction = Fraction(x_final).limit_denominator()
print(f"\nThe exact value for x(1) is {final_fraction.numerator}/{final_fraction.denominator}.")
