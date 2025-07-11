import math

# Define the constants from the derived equation t = sqrt(14 * 350)
factor_from_derivation = 14
height_value = 350

# Calculate the value under the square root
t_squared = factor_from_derivation * height_value

# Calculate the final time
time_to_fall = math.sqrt(t_squared)

# The derived equation for time is t = sqrt(14 * h_0 / g).
# Given h_0 = 350g, the equation becomes t = sqrt(14 * 350g / g).
# The 'g' terms cancel out, leaving t = sqrt(14 * 350).
print(f"The calculation is: t = sqrt({factor_from_derivation} * {height_value})")
print(f"Time to fall: {int(time_to_fall)} seconds")