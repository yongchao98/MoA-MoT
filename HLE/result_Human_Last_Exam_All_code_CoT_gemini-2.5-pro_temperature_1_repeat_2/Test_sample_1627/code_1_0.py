# Define the relationship between the radii.
# We can set r_i to 1 for simplicity, which makes r_o = 2.
r_i = 1
r_o = 2 * r_i

# The formula for the ratio of maximum tangential stress to internal pressure is:
# Ratio = (r_o^2 + r_i^2) / (r_o^2 - r_i^2)

# Calculate the numerator and denominator
numerator = r_o**2 + r_i**2
denominator = r_o**2 - r_i**2

# Print the calculation step-by-step
print(f"The ratio is calculated as (r_o^2 + r_i^2) / (r_o^2 - r_i^2)")
print(f"Given r_o = {r_o} and r_i = {r_i}:")
print(f"Ratio = ({r_o}^2 + {r_i}^2) / ({r_o}^2 - {r_i}^2)")
print(f"Ratio = ({r_o**2} + {r_i**2}) / ({r_o**2} - {r_i**2})")
print(f"Ratio = {numerator} / {denominator}")

# The result is a fraction
print(f"The final ratio is 5/3.")
# To show the decimal value for comparison with answer choice H
print(f"As a decimal, this is approximately: {numerator/denominator:.3f}")