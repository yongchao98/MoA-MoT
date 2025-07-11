import math

# The value r_0 is derived from the analytical solution.
# r_0 = (226 + 49 * sqrt(2)) / 17
# We now calculate the numerical value of r_0.

numerator_c1 = 226
numerator_c2 = 49
denominator = 17

# The equation for r_0
print(f"The equation for the radial distance r_0 is:")
print(f"r_0 = ({numerator_c1} + {numerator_c2} * sqrt(2)) / {denominator}")

# Calculate the value
r_0 = (numerator_c1 + numerator_c2 * math.sqrt(2)) / denominator

# Print the final result
print("\nThe calculated value for r_0 is:")
print(r_0)