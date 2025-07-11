import math

# Given parameters
alpha = 10**16
# R = ln(100/99)
# e^R = 100/99
# e^R - 1 = 1/99
eR_minus_1 = 1/99

# Derived formula for t_0^2
t0_squared = (3 * alpha) / eR_minus_1

# Calculate t_0 (we need the positive value)
t0 = math.sqrt(t0_squared)

# The final equation is t_0 = sqrt( (3 * alpha) / (e^R - 1) )
# Let's print the numbers in the final equation
print("Calculation Steps:")
print(f"alpha = {alpha}")
print(f"R = ln(100/99)")
print(f"e^R - 1 = {eR_minus_1}")
print(f"t_0^2 = (3 * {alpha}) / {eR_minus_1} = {t0_squared}")
print(f"t_0 = sqrt({t0_squared})")
print("\nFinal Value:")
print(f"{t0}")