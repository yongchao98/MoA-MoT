# Define the relationship between outer and inner radius.
# We can set the inner radius to 1 for simplicity, making the outer radius 2.
r_i = 1
r_o = 2 * r_i

# The formula for the ratio of maximum tangential stress (sigma_t_max)
# to internal pressure (p_i) is: (r_o^2 + r_i^2) / (r_o^2 - r_i^2)

# Calculate the numerator and denominator of the ratio
numerator = r_o**2 + r_i**2
denominator = r_o**2 - r_i**2

# Calculate the final ratio
ratio = numerator / denominator

# Print the calculation steps
print("The ratio of maximum tangential stress to internal pressure is calculated as follows:")
print(f"σ_t_max / p_i = (r_o² + r_i²) / (r_o² - r_i²)")
print(f"Given r_o = 2 * r_i, we can use r_i = {r_i} and r_o = {r_o}.")
print(f"Substituting the values into the equation:")
print(f"σ_t_max / p_i = ({r_o}² + {r_i}²) / ({r_o}² - {r_i}²)")
print(f"σ_t_max / p_i = ({r_o**2} + {r_i**2}) / ({r_o**2} - {r_i**2})")
print(f"σ_t_max / p_i = {numerator} / {denominator}")
print(f"The final ratio is {numerator}/{denominator}, which is approximately {ratio:.3f}")
