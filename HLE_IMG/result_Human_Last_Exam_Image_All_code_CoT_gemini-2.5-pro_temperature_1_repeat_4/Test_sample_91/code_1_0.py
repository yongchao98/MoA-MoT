import math

# --- Given parameters ---
L = 40.0  # Length of the drawbridge in meters
y_target = 10.0 # Target vertical height in meters

# --- Step 1: Find sin(theta) at the target height ---
# The relationship is y = L * cos(theta), so cos(theta) = y / L
cos_theta = y_target / L
# Using sin^2(theta) + cos^2(theta) = 1, we find sin(theta)
# Since the bridge is rising, theta is between 0 and pi/2, so sin(theta) is positive.
sin_theta_val = math.sqrt(1 - cos_theta**2)

# --- Step 2: Get the value of dtheta/dt ---
# The problem gives dtheta/dt = -(3*pi/10) / cos(pi/12)
cos_pi_12 = math.cos(math.pi / 12)
dtheta_dt_val = -(3 * math.pi / 10) / cos_pi_12

# --- Step 3: Calculate the vertical velocity (dy/dt) ---
# The formula is dy/dt = -L * sin(theta) * dtheta/dt
dy_dt = -L * sin_theta_val * dtheta_dt_val

# --- Step 4: Output the equation and the final answer ---
print("The equation for the vertical velocity is: dy/dt = -L * sin(theta) * dtheta/dt")
print("\nSubstituting the values into the equation:")
print(f"L = {L}")
print(f"sin(theta) = sqrt(1 - (10/40)^2) = sqrt(15)/4")
print(f"dtheta/dt = -(3*\u03c0/10) / cos(\u03c0/12)")
print("\nSo the full equation with numbers is:")
print("dy/dt = - (40) * (sqrt(15)/4) * (-(3*\u03c0/10) / cos(\u03c0/12))")

print(f"\nThe calculated vertical velocity of the moving edge is: {dy_dt:.4f} m/s")

# For the final answer format
final_answer = dy_dt