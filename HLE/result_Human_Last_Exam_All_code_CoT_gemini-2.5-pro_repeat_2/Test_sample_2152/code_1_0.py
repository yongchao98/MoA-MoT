import math

# Define the constants for the equation
numerator = 3
denominator_factor = 4
angle_rad = math.pi / 30

# Calculate sin(theta)^2
sin_theta_sq = math.sin(angle_rad)**2

# Calculate the final ratio
ratio = numerator / (denominator_factor * sin_theta_sq)

# Output the steps and the final result
print("The ratio of the differential cross-sections is calculated using the formula:")
print("Ratio = Numerator / (Factor * sin(theta)^2)")
print(f"Numerator = {numerator}")
print(f"Factor = {denominator_factor}")
print(f"theta (radians) = pi / 30 = {angle_rad}")
print(f"sin(theta)^2 = sin({angle_rad:.5f})^2 = {sin_theta_sq:.5f}")
print(f"Final Calculation: Ratio = {numerator} / ({denominator_factor} * {sin_theta_sq:.5f})")
print(f"Ratio = {ratio:.3f}")