import math

# Define the given parameters
m = 1.0  # mass of the ring in kg
M = 1.0  # mass of the object in kg
g = 9.8  # acceleration due to gravity in m/s^2
theta_deg = 60.0 # angle in degrees

# Convert angle to radians for use in math functions
theta_rad = math.radians(theta_deg)

# Calculate sine and cosine of the angle
sin_theta = math.sin(theta_rad)
cos_theta = math.cos(theta_rad)

# For clarity in the formula
sin_sq_theta = sin_theta**2
cos_sq_theta = cos_theta**2

# Calculate the numerator of the tension formula
# T_num = m * M * g * sin(theta) * (3*(m+M) - M*sin^2(theta))
T_num = m * M * g * sin_theta * (3 * (m + M) - M * sin_sq_theta)

# Calculate the denominator of the tension formula
# T_den = (m + M*cos^2(theta))^2
T_den = (m + M * cos_sq_theta)**2

# Calculate the tension
Tension = T_num / T_den

# Print the equation with the numerical values
print(f"Calculating the tension T using the formula:")
print(f"T = (m * M * g * sin(θ) * (3*(m+M) - M*sin²(θ))) / (m + M*cos²(θ))²\n")
print(f"Plugging in the values:")
print(f"m = {m} kg")
print(f"M = {M} kg")
print(f"g = {g} m/s²")
print(f"θ = {theta_deg}°\n")

print("Final Equation:")
# Using a formatted string to show the calculation with numbers
equation_str = (
    f"T = ({m} * {M} * {g:.1f} * sin({theta_deg}°) * (3*({m}+{M}) - {M}*sin²({theta_deg}°)))"
    f" / ({m} + {M}*cos²({theta_deg}°))²"
)
print(equation_str)

# Calculate intermediate values for display
part1 = m * M * g * sin_theta
part2 = 3 * (m + M) - M * sin_sq_theta
part3 = (m + M * cos_sq_theta)**2
print(f"T = ({part1:.4f} * {part2:.4f}) / {part3:.4f}")
print(f"T = {T_num:.4f} / {T_den:.4f}\n")

# Print the final result rounded to two decimal places
print(f"The tension in the string is: {Tension:.2f} Newtons")
<<<28.52>>>