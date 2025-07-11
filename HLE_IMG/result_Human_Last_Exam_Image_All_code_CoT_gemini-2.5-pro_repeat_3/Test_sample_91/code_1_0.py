import math

# --- Given information ---
L = 40.0  # Length of the drawbridge in meters
y = 10.0  # Vertical height of the edge of the bridge in meters

# --- Step-by-step solution ---

# The vertical height 'y' of the edge of the drawbridge is related to its length 'L'
# and the angle 'θ' it makes with the vertical by the equation:
# y = L * cos(θ)

# We want to find the vertical velocity, which is the rate of change of y with respect to time, dy/dt.
# We can find this by differentiating the equation with respect to time 't' using the chain rule:
# dy/dt = d/dt [L * cos(θ)]
# dy/dt = -L * sin(θ) * dθ/dt

# --- Calculations ---

# Step 1: Find the values of cos(θ) and sin(θ) when y = 10 m.
# From y = L * cos(θ), we get cos(θ) = y / L.
cos_theta = y / L

# Using the identity sin²(θ) + cos²(θ) = 1, we find sin(θ).
# Since the bridge is being raised from the horizontal, θ is between 0 and π/2, so sin(θ) is positive.
sin_theta = math.sqrt(1 - cos_theta**2)

# Step 2: Use the given rate of change of the angle, dθ/dt.
# The problem states dθ/dt = - (3 * π / 10) / cos(π / 12)
# Let's calculate the numerical values for the components of dθ/dt.
pi = math.pi
cos_pi_over_12 = math.cos(pi / 12)
dtheta_dt = -(3 * pi / 10) / cos_pi_over_12

# Step 3: Substitute all the calculated values into the equation for dy/dt.
# dy/dt = -L * sin(θ) * dθ/dt
vertical_velocity = -L * sin_theta * dtheta_dt

# --- Output the results ---

print("Problem: Find the vertical velocity of the drawbridge's edge.")
print("\n1. The equation for the vertical velocity (dy/dt) is derived from y = L * cos(θ):")
print("   dy/dt = -L * sin(θ) * dθ/dt")

print("\n2. We evaluate each term in the equation for when the height y = 10 m:")
print(f"   - L (length of bridge) = {L} m")

print(f"\n   - To find sin(θ):")
print(f"     cos(θ) = y / L = {y} / {L} = {cos_theta}")
print(f"     sin(θ) = sqrt(1 - cos²(θ)) = sqrt(1 - {cos_theta}**2) = {sin_theta:.6f}")

print(f"\n   - The given rate of change of the angle is dθ/dt = -(3 * π / 10) / cos(π / 12):")
print(f"     3 * π / 10 = {3 * pi / 10:.6f}")
print(f"     cos(π / 12) = {cos_pi_over_12:.6f}")
print(f"     dθ/dt = -({3 * pi / 10:.6f}) / {cos_pi_over_12:.6f} = {dtheta_dt:.6f} rad/s")

print("\n3. Substitute these numbers into the velocity equation:")
print(f"   dy/dt = -({L}) * ({sin_theta:.6f}) * ({dtheta_dt:.6f})")

print(f"\n4. Final Calculation:")
print(f"   dy/dt = {vertical_velocity:.4f} m/s")