import math

# Define the exponents from the problem statement
alpha_exp = 10000
x0_exp = -5000000
y0_exp = -5000000

# The problem simplifies to finding T using the first-order solvability condition.
# The derived equation for T is T ≈ α / (x₀ + y₀), as x₀ and y₀ are very small.
# We can calculate T by computing the exponents.

# Calculate the final exponent for T
# T ≈ (10^alpha_exp) / (10^x0_exp + 10^y0_exp)
# T ≈ (10^alpha_exp) / (2 * 10^x0_exp)
# T ≈ 0.5 * 10^(alpha_exp - x0_exp)
final_exponent = alpha_exp - x0_exp

# Create string representations for printing the equation
alpha_str = f"10^({alpha_exp})"
x0_str = f"10^({x0_exp})"
y0_str = f"10^({y0_exp})"
denominator_str = f"2 * 10^({x0_exp})"
final_T_str = f"0.5 * 10^({final_exponent})"

print("The simplified equation for T is derived from the solvability condition:")
print(f"T = α / (x₀ + y₀)")
print("\nSubstituting the given values:")
print(f"T = {alpha_str} / ({x0_str} + {y0_str})")
print(f"T = {alpha_str} / ({denominator_str})")
print("\nThe final calculated value for T is:")
print(f"T = {final_T_str}")