import math

# We are given the condition x(0) > 1. Let's choose a sample value.
# x_0 = 8 is a convenient choice because 8^(2/3) is an integer.
x0 = 8.0

# The condition for blow-up is y(0) < sqrt(2*x(0) + 1 - 3*x(0)^(2/3)).
# We will now calculate the right-hand side of this inequality.

# Calculate the value under the square root
val_2_x0 = 2 * x0
val_x0_pow_2_3 = x0**(2.0/3.0)
val_3_x0_pow_2_3 = 3 * val_x0_pow_2_3
rhs_value_sq = val_2_x0 + 1 - val_3_x0_pow_2_3
final_threshold = math.sqrt(rhs_value_sq)

# Print the step-by-step calculation
print(f"To find the values of y(0) for which the solution blows up, we use the inequality:")
print(f"y(0) < sqrt(2*x(0) + 1 - 3*x(0)^(2/3))")
print(f"\nGiven x(0) = {x0:.0f}, we substitute this value into the equation:")
print(f"y(0) < sqrt(2 * {x0:.0f} + 1 - 3 * ({x0:.0f})^(2/3))")
# Show the numbers in the equation
print(f"y(0) < sqrt({val_2_x0:.0f} + 1 - 3 * {val_x0_pow_2_3:.0f})")
print(f"y(0) < sqrt({val_2_x0:.0f} + 1 - {val_3_x0_pow_2_3:.0f})")
print(f"y(0) < sqrt({rhs_value_sq:.0f})")
print(f"\nTherefore, the solution blows up for y(0) < {final_threshold}")
