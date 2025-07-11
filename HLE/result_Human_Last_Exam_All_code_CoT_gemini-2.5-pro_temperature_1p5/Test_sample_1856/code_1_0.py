import math

# Given values
beta_A = 0.95

# Calculate the Lorentz factor gamma for particle A
gamma_A = 1 / math.sqrt(1 - beta_A**2)

# Value of sqrt(2)
sqrt2 = math.sqrt(2)

# Calculate the tangent of the angle theta_C using the derived formula
# tan(theta_C) = 1 / (gamma_A * (1 + sqrt(2) * beta_A))
denominator = gamma_A * (1 + sqrt2 * beta_A)
tan_theta = 1 / denominator

# Calculate the angle in radians
angle_rad = math.atan(tan_theta)

# Convert the angle to degrees
angle_deg = math.degrees(angle_rad)

# Print the formula and the values used
print("The formula for the tangent of the angle is:")
print("tan(theta) = 1 / (gamma_A * (1 + sqrt(2) * beta_A))")
print("\nPlugging in the numerical values:")
print(f"beta_A = {beta_A}")
print(f"gamma_A = 1 / sqrt(1 - {beta_A}^2) = {gamma_A:.4f}")
print(f"sqrt(2) = {sqrt2:.4f}")
print("\nSo the equation becomes:")
print(f"tan(theta) = 1 / ({gamma_A:.4f} * (1 + {sqrt2:.4f} * {beta_A}))")
print(f"tan(theta) = {tan_theta:.4f}")

# Print the final result
print(f"\nThe angle between the direction of C and that of its parent particle A is: {angle_deg:.3f} degrees.")
