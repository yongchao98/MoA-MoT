import math

# Given the velocity of particle A
beta_A = 0.95

# Calculate the Lorentz factor, gamma
gamma_A = 1 / math.sqrt(1 - beta_A**2)

# The angle theta_C is given by arctan(1 / (gamma_A * (1 + sqrt(2) * beta_A)))
# Let's calculate the terms in the equation.
term_in_parentheses = 1 + math.sqrt(2) * beta_A
denominator = gamma_A * term_in_parentheses
value_inside_arctan = 1 / denominator

# Calculate the angle in radians
theta_rad = math.atan(value_inside_arctan)

# Convert the angle to degrees
theta_deg = math.degrees(theta_rad)

# Print the final equation with the calculated numbers
print("The final equation for the angle theta_C is:")
print(f"theta_C = arctan(1 / (gamma_A * (1 + sqrt(2) * beta_A)))")
print("\nPlugging in the numbers:")
print(f"beta_A = {beta_A}")
print(f"gamma_A = 1 / sqrt(1 - {beta_A}^2) = {gamma_A:.4f}")
print(f"theta_C = arctan(1 / ({gamma_A:.4f} * (1 + sqrt(2) * {beta_A})))")
print(f"theta_C = arctan(1 / ({gamma_A:.4f} * {term_in_parentheses:.4f}))")
print(f"theta_C = arctan(1 / {denominator:.4f})")
print(f"theta_C = arctan({value_inside_arctan:.4f})")
print(f"\nThe final angle is {theta_deg:.3f} degrees.")
