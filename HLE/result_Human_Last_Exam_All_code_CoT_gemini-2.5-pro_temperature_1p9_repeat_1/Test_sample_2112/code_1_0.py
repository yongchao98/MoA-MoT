import math

# These are the numbers from the derived final equation for r_0.
num_const = 226
sqrt_coeff = 49
sqrt_val = 2
denominator = 17

# Calculate the value of r_0 using the derived expression.
r_0 = (num_const + sqrt_coeff * math.sqrt(sqrt_val)) / denominator

# Print the final equation with each of its constituent numbers.
print(f"The final equation for r_0 is:")
print(f"r_0 = ({num_const} + {sqrt_coeff} * sqrt({sqrt_val})) / {denominator}")

# Print the final numerical answer.
print("\nThe calculated value of the radial distance r_0 is:")
print(r_0)