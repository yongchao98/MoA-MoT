import math

# Define the given physical constants
m = 0.1   # mass of the block in kg (100 g)
M = 10.0  # mass of the wedge in kg
h = 2.0   # height of the wedge in m
theta_deg = 30.0 # angle of the wedge in degrees
mu = 0.5  # coefficient of friction
g = 10.0  # acceleration due to gravity in m/s^2

# Convert angle from degrees to radians for use in math functions
theta_rad = math.radians(theta_deg)

# --- Calculation ---
# Step 1: Calculate sin(theta) and cos(theta)
sin_theta = math.sin(theta_rad)
cos_theta = math.cos(theta_rad)

# Step 2: Calculate the relative acceleration (a_rel)
# a_rel = [g * (M + m) * (sin(theta) - mu * cos(theta))] / [M + m * sin(theta) * (sin(theta) - mu * cos(theta))]

# Calculate the term (sin(theta) - mu * cos(theta)) which appears twice
common_term = sin_theta - mu * cos_theta

# Calculate numerator of a_rel
a_rel_num = g * (M + m) * common_term

# Calculate denominator of a_rel
a_rel_den = M + m * sin_theta * common_term

# Calculate a_rel
a_rel = a_rel_num / a_rel_den

print("--- Step 1: Calculate the relative acceleration (a_rel) of the block ---")
print("Using formula: a_rel = [g*(M+m)*(sin(θ)-μ*cos(θ))] / [M+m*sin(θ)*(sin(θ)-μ*cos(θ))]")
print(f"Numerator = g * (M + m) * (sin({theta_deg}) - μ * cos({theta_deg}))")
print(f"Numerator = {g} * ({M} + {m}) * ({sin_theta:.4f} - {mu} * {cos_theta:.4f}) = {a_rel_num:.4f} N")
print(f"Denominator = M + m * sin({theta_deg}) * (sin({theta_deg}) - μ * cos({theta_deg}))")
print(f"Denominator = {M} + {m} * {sin_theta:.4f} * ({sin_theta:.4f} - {mu} * {cos_theta:.4f}) = {a_rel_den:.4f} kg")
print(f"a_rel = {a_rel_num:.4f} / {a_rel_den:.4f} = {a_rel:.4f} m/s²\n")

# Step 3: Calculate the length of the slope (L)
# L = h / sin(theta)
L = h / sin_theta

print("--- Step 2: Calculate the distance (L) the block slides down ---")
print(f"L = h / sin(θ) = {h} / sin({theta_deg}°) = {h} / {sin_theta:.4f} = {L:.4f} m\n")

# Step 4: Calculate the time (t)
# t = sqrt(2 * L / a_rel)
if a_rel > 0:
    t = math.sqrt(2 * L / a_rel)
    print("--- Step 3: Calculate the time (t) to slide down ---")
    print("t = sqrt(2 * L / a_rel)")
    print(f"t = sqrt(2 * {L:.4f} / {a_rel:.4f})")
    print(f"The time it takes for the block to slide to the bottom is: {t:.4f} seconds.")
else:
    print("The acceleration is not positive, meaning the block will not slide down due to high friction.")

# Final answer in the required format
print(f"\n<<<{t:.4f}>>>")
