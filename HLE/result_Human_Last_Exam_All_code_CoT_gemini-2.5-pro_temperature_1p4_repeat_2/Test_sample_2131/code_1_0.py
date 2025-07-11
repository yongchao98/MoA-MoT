import math

# This script calculates the membrane's deflection at x = 0, which is y(0).
# Based on the analytical solution of the differential equation,
# the value of y(0) is given by the formula:
# y(0) = (1/3) * (9/4)^(8/5)

# Define the numbers in the final equation
c = 1/3
base_num = 9
base_den = 4
exp_num = 8
exp_den = 5

base = base_num / base_den
exponent = exp_num / exp_den

# Calculate the result
y_0 = c * (base ** exponent)

# Output the equation and the final result
print("The deflection y(0) is calculated using the formula:")
print(f"y(0) = ({c:.1f}) * ({base_num}/{base_den})**({exp_num}/{exp_den})")
print("\nCalculating the values:")
print(f"Base of the exponent = {base}")
print(f"Exponent = {exponent}")
print(f"y(0) = {c:.6f} * {base:.2f}**{exponent:.1f}")
print(f"y(0) = {c:.6f} * {base ** exponent:.6f}")
print(f"\nThe membrane's deflection at x = 0 is y(0) = {y_0}")
