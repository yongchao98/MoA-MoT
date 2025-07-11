import math

# Define the given physical constants
m = 1.0  # Mass of the ring in kg
M = 1.0  # Mass of the object in kg
g = 9.8  # Acceleration due to gravity in m/s^2
theta_deg = 60.0 # Angle of the string in degrees

# Convert the angle from degrees to radians for use in trigonometric functions
theta_rad = math.radians(theta_deg)

# For clarity in the final equation, we calculate sin(theta) and cos(theta) first
sin_theta = math.sin(theta_rad)
cos_theta = math.cos(theta_rad)

# The tension 'T' can be derived from the principles of energy and momentum conservation.
# The resulting formula for tension is:
# T = (m * M * g * sin(theta) * (3*m + 2*M + M*cos(theta)^2)) / (m + M*cos(theta)^2)^2

print("This script calculates the tension using the derived formula:")
print("T = (m*M*g*sin(theta) * (3*m + 2*M + M*(cos(theta))^2)) / (m + M*(cos(theta))^2)^2\n")
print("--- Calculation Breakdown ---")

# Calculate the value of each term in the formula
# Numerator parts
numerator_part1 = m * M * g * sin_theta
numerator_part2 = 3 * m + 2 * M + M * cos_theta**2
numerator = numerator_part1 * numerator_part2

# Denominator part
denominator_base = m + M * cos_theta**2
denominator = denominator_base**2

# Final Tension
tension = numerator / denominator

# Print the values of each number used in the final equation
print("Plugging in the given values:")
print(f"m = {m} kg")
print(f"M = {M} kg")
print(f"g = {g} m/s^2")
print(f"theta = {theta_deg} degrees")
print(f"sin(theta) = {sin_theta:.4f}")
print(f"cos(theta)^2 = {cos_theta**2:.4f}\n")

print(f"Numerator = ({m} * {M} * {g} * {sin_theta:.4f}) * (3*{m} + 2*{M} + {M}*{cos_theta**2:.4f})")
print(f"Numerator = ({numerator_part1:.4f}) * ({numerator_part2:.4f}) = {numerator:.4f}\n")

print(f"Denominator = ({m} + {M} * {cos_theta**2:.4f})^2")
print(f"Denominator = ({denominator_base:.4f})^2 = {denominator:.4f}\n")

print(f"Tension (T) = {numerator:.4f} / {denominator:.4f}\n")

print(f"The calculated tension is {tension:.4f} Newtons.")
print("When rounded to two decimal places, the tension is:")
print(f"{tension:.2f}")

<<<28.52>>>