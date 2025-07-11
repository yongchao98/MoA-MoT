import math

# --- Given parameters ---
m = 0.1  # kg (mass of the block, 100g)
M = 10.0 # kg (mass of the wedge)
theta_deg = 30.0 # degrees (angle of the wedge)
h = 2.0  # m (height of the wedge)
mu = 0.5 # coefficient of kinetic friction
g = 10.0 # m/s^2 (acceleration due to gravity)

# --- Step 1: Convert angle to radians and calculate slope length ---
theta_rad = math.radians(theta_deg)
L = h / math.sin(theta_rad)

# --- Step 2: Solve for accelerations using the non-inertial frame method ---
# We solve the system of equations derived from Newton's laws for the block and wedge.
# Let's define some shorthand for trigonometric functions
S = math.sin(theta_rad)
C = math.cos(theta_rad)

# Equation for Normal Force (N)
# Derived from summing forces perpendicular to the slope in the wedge's frame
# and forces on the wedge in the ground frame.
# N = (m*g*C) / (1 + (m/M)*(S**2 + mu*S*C))
N_numerator = m * g * C
N_denominator = 1 + (m / M) * (S**2 + mu * S * C)
N = N_numerator / N_denominator

# Equation for Wedge Acceleration (Ax)
# Derived from summing horizontal forces on the wedge.
# Ax = N/M * (S + mu*C)
Ax = (N / M) * (S + mu * C)

# Equation for Block's Relative Acceleration (a_rel)
# Derived from summing forces parallel to the slope in the wedge's frame.
# a_rel = g*S - (mu/m)*N - Ax*C
a_rel = g * S - (mu / m) * N - Ax * C

# --- Step 3: Calculate the time (t) using kinematics ---
# L = v0*t + 0.5*a_rel*t^2, where v0 (initial relative velocity) is 0.
# t = sqrt(2 * L / a_rel)
t_squared = (2 * L) / a_rel
t = math.sqrt(t_squared)

# --- Step 4: Output the results ---
print(f"Problem Parameters:")
print(f"Block mass (m) = {m} kg")
print(f"Wedge mass (M) = {M} kg")
print(f"Wedge angle (theta) = {theta_deg} degrees")
print(f"Wedge height (h) = {h} m")
print(f"Friction coefficient (mu) = {mu}")
print(f"Gravity (g) = {g} m/s^2\n")

print("--- Calculation ---")
print(f"The distance the block slides is L = h / sin(theta) = {L:.4f} m.")
print(f"The relative acceleration of the block down the slope is a_rel = {a_rel:.4f} m/s^2.")
print("\nThe time 't' is found using the kinematic equation: t = sqrt(2 * L / a_rel)")
print(f"Substituting the values: t = sqrt(2 * {L:.4f} / {a_rel:.4f})")
print(f"The final calculated time is: {t:.4f} seconds.")

# Final answer in the required format
final_answer = t
print(f"\n<<<{final_answer}>>>")