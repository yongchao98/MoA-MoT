import math

# Define the velocity of particle A
beta_A = 0.95

# Calculate the Lorentz factor gamma_A
gamma_A = 1 / math.sqrt(1 - beta_A**2)

# Calculate the tangent of the angle theta_C in the lab frame
# The formula is tan(theta_C) = (1/sqrt(2)) / (gamma_A * (1/sqrt(2) + beta_A))
numerator = 1 / math.sqrt(2)
denominator = gamma_A * (1 / math.sqrt(2) + beta_A)
tan_theta_C = numerator / denominator

# Calculate the angle in radians and then convert to degrees
theta_C_rad = math.atan(tan_theta_C)
theta_C_deg = math.degrees(theta_C_rad)

# Print the breakdown of the final calculation as requested
print("The equation for the tangent of the angle is: tan(theta_C) = (1/sqrt(2)) / (gamma_A * (1/sqrt(2) + beta_A))")
print("\nSubstituting the values:")
print(f"beta_A = {beta_A}")
print(f"gamma_A = 1 / sqrt(1 - {beta_A}^2) = {gamma_A:.5f}")
print(f"sqrt(2) = {math.sqrt(2):.5f}")

print("\nEvaluating the equation:")
print(f"tan(theta_C) = ({numerator:.5f}) / ({gamma_A:.5f} * ({numerator:.5f} + {beta_A}))")
print(f"tan(theta_C) = {numerator:.5f} / ({denominator:.5f})")
print(f"tan(theta_C) = {tan_theta_C:.5f}")

# Print the final result
print(f"\nThe angle theta_C is arctan({tan_theta_C:.5f}) which is {theta_C_deg:.3f} degrees.")
print(f"\nFinal Answer (rounded to three decimal places):")
print(f"{theta_C_deg:.3f}")