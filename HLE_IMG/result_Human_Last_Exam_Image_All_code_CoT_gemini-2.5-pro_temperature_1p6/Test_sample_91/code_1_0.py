import math

# Given parameters from the problem
L = 40.0  # Length of the drawbridge in meters
h_target = 10.0   # Vertical height of the edge at the moment of interest in meters

# --- Step-by-step calculation ---

# 1. Find cos(theta) at the specified height
# The relationship is h = L * cos(theta)
cos_theta = h_target / L

# 2. Find sin(theta) from cos(theta)
# Using sin^2(theta) + cos^2(theta) = 1
# sin(theta) must be positive as the bridge is rising (0 < theta < pi/2)
sin_theta = math.sqrt(1 - cos_theta**2)

# 3. Calculate the given rate of change of theta, d(theta)/dt
# d(theta)/dt = -3*pi / (10 * cos(pi/12))
dtheta_dt = -3 * math.pi / (10 * math.cos(math.pi / 12))

# 4. Calculate the vertical velocity, dh/dt
# The formula is dh/dt = -L * sin(theta) * d(theta)/dt
dh_dt = -L * sin_theta * dtheta_dt

# --- Output the results ---
print("This script calculates the vertical velocity of the drawbridge edge.")
print("\nThe final velocity is calculated using the equation: dh/dt = -L * sin(theta) * d(theta)/dt\n")
print("Here are the values for each component of the equation:")
print(f"L (Length of bridge) = {L}")
print(f"sin(theta) [when height is {h_target}m] = {sin_theta}")
print(f"d(theta)/dt (Rate of change of angle) = {dtheta_dt} rad/s")

# Final result
print(f"\nFinal calculated vertical velocity (dh/dt) = {dh_dt} m/s")