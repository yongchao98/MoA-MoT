#
# Plan:
# 1. Define the relationship between the outer and inner radii. The problem states r_o = 2 * r_i.
# 2. For simplicity, we can assume r_i = 1, which makes r_o = 2. The ratio is independent of the actual radius value.
# 3. The formula for the ratio of maximum tangential stress to internal pressure is (r_o^2 + r_i^2) / (r_o^2 - r_i^2).
# 4. Calculate the numerator and denominator using the assumed radii values.
# 5. Compute the final ratio.
# 6. Print the full equation with the calculated values.
#

# Let the inner radius (r_i) be 1 unit and the outer radius (r_o) be 2 units.
r_i = 1
r_o = 2

# The formula for the ratio of maximum tangential stress (sigma_t_max)
# to internal pressure (P_i) is: (r_o^2 + r_i^2) / (r_o^2 - r_i^2)

# Calculate the numerator part of the equation
numerator_val = r_o**2 + r_i**2

# Calculate the denominator part of the equation
denominator_val = r_o**2 - r_i**2

# Calculate the final ratio
ratio = numerator_val / denominator_val

# Print the step-by-step calculation
print(f"The equation for the ratio is: (r_o^2 + r_i^2) / (r_o^2 - r_i^2)")
print(f"Substituting r_o = {r_o} and r_i = {r_i}:")
print(f"Ratio = ({r_o}^2 + {r_i}^2) / ({r_o}^2 - {r_i}^2)")
print(f"Ratio = ({r_o**2} + {r_i**2}) / ({r_o**2} - {r_i**2})")
print(f"Ratio = {numerator_val} / {denominator_val}")
print(f"The numerical value of the ratio is {ratio}")
print(f"As a fraction, the ratio is {numerator_val}/{denominator_val}.")
