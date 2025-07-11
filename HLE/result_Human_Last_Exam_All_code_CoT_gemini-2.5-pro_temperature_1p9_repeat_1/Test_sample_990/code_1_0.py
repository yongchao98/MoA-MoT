import math

# The derived formula for the time of fall T is T = sqrt(14 * h / g).
# h is the height and g is the acceleration due to gravity.

# The height is given as h = 350g meters.
h_factor = 350

# Substitute h = 350g into the formula: T = sqrt(14 * (350 * g) / g)
# The variable g cancels out.
final_value_under_sqrt = 14 * h_factor
result = math.sqrt(final_value_under_sqrt)

print("The final equation for the time T is T = sqrt(14 * h / g).")
print(f"Given the height h = {h_factor}*g meters, we substitute this into the equation.")
print(f"T = sqrt((14 * ({h_factor} * g)) / g)")
print(f"The variable 'g' cancels out, leaving:")
print(f"T = sqrt(14 * {h_factor})")
print(f"T = sqrt({final_value_under_sqrt})")
print(f"T = {result} seconds")