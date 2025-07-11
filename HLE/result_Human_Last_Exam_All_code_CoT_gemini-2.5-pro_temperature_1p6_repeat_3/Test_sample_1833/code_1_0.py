import math

# --- Given values ---
# Mass of the ring and the object (m = M)
m = 1.0  # kg
# Acceleration due to gravity
g = 9.8  # m/s^2
# Angle in degrees
theta_deg = 60.0 # degrees

# --- Calculations ---
# Convert angle to radians for use in math functions
theta_rad = math.radians(theta_deg)

# Calculate sine and cosine of the angle
sin_theta = math.sin(theta_rad)
cos_theta = math.cos(theta_rad)

# The formula for tension T when m = M is:
# T = (m * g * sin(theta) * (5 + cos(theta)^2)) / (1 + cos(theta)^2)^2

# Numerator and denominator terms
numerator_val = m * g * sin_theta * (5 + cos_theta**2)
denominator_val = (1 + cos_theta**2)**2

# Calculate the tension
tension = numerator_val / denominator_val

# --- Output the result ---
# The problem asks to show the numbers in the final equation.
print("Calculating the tension in the string.")
print("The formula for tension T, when m=M, is:")
print("T = (m * g * sin(θ) * (5 + cos(θ)²)) / (1 + cos(θ)² )²\n")

print("Plugging in the values:")
print(f"m = {m} kg")
print(f"g = {g} m/s²")
print(f"θ = {theta_deg}°\n")

# Show the equation with the numerical values
print("Equation with numbers:")
print(f"T = ({m} * {g} * sin({theta_deg}°) * (5 + cos({theta_deg}°)²)) / (1 + cos({theta_deg}°)² )²")
print(f"T = ({m} * {g:.1f} * {sin_theta:.3f} * (5 + {cos_theta:.3f}²)) / (1 + {cos_theta:.3f}² )²")
print(f"T = ({m * g * sin_theta:.3f} * {5 + cos_theta**2:.3f}) / ({1 + cos_theta**2:.3f})²")
print(f"T = {numerator_val:.3f} / {denominator_val:.3f}")
print(f"T = {tension:.2f} N\n")

# The final answer as required by the format
print("The tension in the string is calculated below, rounded to two decimal places.")
print(f"{tension:.2f}")