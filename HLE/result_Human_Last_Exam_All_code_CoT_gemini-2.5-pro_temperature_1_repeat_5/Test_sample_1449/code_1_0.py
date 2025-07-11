import math

# --- Given Parameters ---
m = 100 / 1000  # Mass of the block in kg
M = 10.0        # Mass of the wedge in kg
theta_deg = 30.0# Angle of the wedge in degrees
h = 2.0         # Height of the wedge in m
mu = 0.5        # Coefficient of kinetic friction
g = 10.0        # Acceleration due to gravity in m/s^2

# Note: The initial velocity v0 = 1 m/s is not needed to find the time taken
# for the block to slide down relative to the wedge, as accelerations are
# independent of the initial velocity of the system.

# --- Calculations ---

# Convert angle to radians for use in trigonometric functions
theta_rad = math.radians(theta_deg)
sin_t = math.sin(theta_rad)
cos_t = math.cos(theta_rad)

# To find the time, we first need the acceleration of the block relative
# to the wedge (a_rel). We can derive this by analyzing the forces on both
# the block and the wedge.

# Let 'A' be the magnitude of the wedge's horizontal acceleration.
# The formula for 'A', derived from Newton's laws, is:
# A = (m*g*(sin(t)*cos(t) + mu*cos(t)**2)) / (M - m*sin(t)**2 - m*mu*sin(t)*cos(t))
numerator_A = m * g * (sin_t * cos_t + mu * cos_t**2)
denominator_A = M - m * sin_t**2 - m * mu * sin_t * cos_t
A = numerator_A / denominator_A

# The acceleration of the block relative to the wedge (a_rel) is given by:
# a_rel = g*(sin(t) - mu*cos(t)) + A*(cos(t) - mu*sin(t))
a_rel = g * (sin_t - mu * cos_t) + A * (cos_t - mu * sin_t)

# The distance the block slides along the incline is L.
L = h / sin_t

# The time 't' is found using the kinematic equation L = 0.5 * a_rel * t^2.
# We must check that a_rel is positive, meaning the block will slide down.
if a_rel <= 0:
    print("The block does not slide down; friction is too strong.")
    time = float('inf')
else:
    time = math.sqrt(2 * L / a_rel)

# --- Output the Final Answer ---

print("The time 't' for the block to slide down the wedge is found using the equation:")
print("t = sqrt(2 * L / a_rel)")
print("\nFirst, we calculate the values for L and a_rel:")
print(f"The distance along the incline, L = {h} / sin({theta_deg}) = {L:.4f} m")
print(f"The relative acceleration of the block, a_rel = {a_rel:.4f} m/s^2")
print("\nSubstituting these values into the equation:")
print(f"t = sqrt(2 * {L:.4f} / {a_rel:.4f})")
print(f"The calculated time is: {time:.4f} s")

# The final answer in the required format
# print(f'<<<{time:.4f}>>>')