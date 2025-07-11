import math

# Step 1: Define physical constants and geometry
m = 0.1  # kg, mass of the block
M = 10.0 # kg, mass of the wedge
g = 10.0 # m/s^2, acceleration due to gravity
h = 2.0  # m, initial height of the block
theta_deg = 30.0 # degrees, angle of the wedge
mu = 0.5 # coefficient of friction

# The initial velocity v0 = 1 m/s is not needed for this calculation, 
# as the accelerations are independent of the initial velocity of the system.

# Convert angle to radians for trigonometric functions
theta_rad = math.radians(theta_deg)

# Calculate the distance 'd' the block slides along the incline
d = h / math.sin(theta_rad)

# Step 2: Determine the accelerations
# First, calculate the horizontal acceleration of the wedge, A_x.
# This formula is derived from the equations of motion for the block and wedge.
s = math.sin(theta_rad)
c = math.cos(theta_rad)
s_sq = s**2
c_sq = c**2
s_c = s * c

# Numerator of the expression for A_x
A_x_num = g * m * (mu * c_sq - s_c)

# Denominator of the expression for A_x
A_x_den = M - m * s_sq + m * mu * s_c

# Calculate A_x
A_x = A_x_num / A_x_den

# Next, calculate the acceleration of the block relative to the wedge, a_rel.
# This is derived from the kinematic relationship between A_x and a_rel.
a_rel = -A_x * (M + m) / (m * c)

# Step 3: Calculate the time 't'
# The block starts from rest relative to the wedge, so we use d = (1/2) * a_rel * t^2.
# The final equation for time is t = sqrt(2 * d / a_rel).
t_sq = 2 * d / a_rel
t = math.sqrt(t_sq)

# Step 4: Output the results as requested
print(f"The calculation is based on the following values:")
print(f"Distance to slide, d = {d:.4f} m")
print(f"Relative acceleration, a_rel = {a_rel:.4f} m/s^2\n")

print(f"The final equation for time is: t = sqrt(2 * d / a_rel)")
print(f"Plugging in the numbers: t = sqrt(2 * {d:.4f} / {a_rel:.4f})")
print(f"The calculated time is t = {t:.4f} s.")
print(f"<<<{t}>>>")