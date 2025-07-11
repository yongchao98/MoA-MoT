import math

# --- Given values ---
m = 1.0  # mass of the ring in kg
M = 1.0  # mass of the object in kg
theta_deg = 60.0  # angle in degrees
g = 9.8  # acceleration due to gravity in m/s^2

# Convert angle from degrees to radians for math functions
theta_rad = math.radians(theta_deg)

# Calculate sine and cosine of the angle
sin_theta = math.sin(theta_rad)
cos_theta = math.cos(theta_rad)

# --- Calculate tension using the derived formula ---
# Numerator of the formula: m * M * g * sin(θ) * (3m + M*(2 + cos(θ)²))
numerator_val = m * M * g * sin_theta * (3 * m + M * (2 + cos_theta**2))

# Denominator of the formula: (m + M*cos(θ)²)²
denominator_val = (m + M * cos_theta**2)**2

# Tension T
tension = numerator_val / denominator_val

# --- Output the results ---
print("Calculation of Tension in the String")
print("-" * 35)
print("The formula for tension T is:")
print("T = (m * M * g * sin(θ) * (3*m + M*(2 + cos(θ)²))) / (m + M*cos(θ)²)²\n")

print("Plugging in the given values:")
print(f"Mass of ring (m) = {m} kg")
print(f"Mass of object (M) = {M} kg")
print(f"Angle (θ) = {theta_deg}°")
print(f"Gravity (g) = {g} m/s²\n")

print("The equation with numerical values is:")
# The format below shows each number in the equation.
print(f"T = ({m} * {M} * {g} * {sin_theta:.4f} * (3*{m} + {M}*(2 + {cos_theta:.4f}²))) / ({m} + {M}*{cos_theta:.4f}²)²")

print(f"\nCalculated numerator = {numerator_val:.4f}")
print(f"Calculated denominator = {denominator_val:.4f}")

print(f"\nTension (T) = {numerator_val:.4f} / {denominator_val:.4f}")
print(f"Tension (T) = {tension:.2f} Newtons")

print("\n<<<" + f"{tension:.2f}" + ">>>")