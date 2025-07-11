import math

# Step 1: Define all given constants and calculate key time points and positions from the man's journey.
# These values are taken directly from the problem statement.
u = 1.5  # m/s, man's initial and later constant speed
t0 = 0.0
t1 = 4.0  # s
alpha_obtuse_deg = 130.0  # degrees
a1 = -0.15  # m/s^2, man's deceleration
t3 = 15.0  # s
a2 = 0.25  # m/s^2, man's acceleration
gamma_deg = 40.0  # degrees
t5 = 23.0  # s
t6 = 40.0  # s

# There is a contradiction in the problem statement for the first segment.
# The bird is said to have "northward" motion (positive y-component), but the projection angle of 130 degrees
# with the North axis implies a negative y-component (cos(130) < 0).
# We assume the intended angle between the paths is the acute angle, 180 - 130 = 50 degrees.
# This makes the "northward" and "eastward" descriptions physically consistent.
alpha_deg = 180.0 - alpha_obtuse_deg # 50 degrees
alpha_rad = math.radians(alpha_deg)
gamma_rad = math.radians(gamma_deg)

# Determine key time points and positions for the man based on the original plan.
# t2 is when the man comes to a stop after decelerating.
t2 = t1 + (-u / a1)  # v_final = v_initial + a*t => 0 = u + a1*(t2-t1)
# t4 is when the man accelerates from rest back to speed u.
t4 = t3 + (u / a2) # u = 0 + a2*(t4-t3)

# Man's y-positions (northward displacement) at these key times.
y_man_1 = u * (t1 - t0)
y_man_2 = y_man_1 + u * (t2 - t1) + 0.5 * a1 * (t2 - t1)**2
y_man_3 = y_man_2 # Man stays still between t2 and t3.
y_man_4 = y_man_3 + 0.5 * a2 * (t4 - t3)**2 # Starts from rest at t3.

# Step 2: Solve for the bird's velocity components relative to its total speed, v.
# Let r_x, r_y, r_z be the ratios of the bird's initial velocity components to its speed v.
# From the original plan's rendezvous conditions, we can form a system of equations.
# This leads to a quadratic equation for r_z.
K1 = math.sin(alpha_rad) * math.sin(gamma_rad)
K2 = 1.5 * math.cos(gamma_rad) - 2.5 * math.sin(gamma_rad)
A_quad = K1**2 + math.cos(gamma_rad)**2
B_quad = 2 * math.cos(gamma_rad) * K2
C_quad = K2**2 - K1**2

# Solve the quadratic equation for r_z.
discriminant = B_quad**2 - 4 * A_quad * C_quad
r_z1 = (-B_quad + math.sqrt(discriminant)) / (2 * A_quad)
r_z2 = (-B_quad - math.sqrt(discriminant)) / (2 * A_quad)

# A valid solution for r_z must be positive ("upward" motion) and satisfy the pre-squaring equation.
r_z = 0
if r_z1 > 0:
    lhs = K1 * math.sqrt(1 - r_z1**2)
    rhs = r_z1 * math.cos(gamma_rad) + K2
    if abs(lhs - rhs) < 1e-9:
        r_z = r_z1
if r_z == 0 and r_z2 > 0:
    lhs = K1 * math.sqrt(1 - r_z2**2)
    rhs = r_z2 * math.cos(gamma_rad) + K2
    if abs(lhs - rhs) < 1e-9:
        r_z = r_z2

# Calculate the other velocity ratios.
r_x = math.sin(alpha_rad) * math.sqrt(1 - r_z**2)
r_y = r_x * math.cot(alpha_rad)

# Step 3: Solve for the bird's constant speed, v.
delta_t_orig = (4 * r_z + 6) / math.sin(gamma_rad)
v_numerator = y_man_4 + u * delta_t_orig
v_denominator = 4 * r_y + 1
v = v_numerator / v_denominator

# Step 4: Calculate the positions of the man and the bird at t5 = 23s.
y_man_5 = y_man_4 + u * (t5 - t4)

# To find the bird's position at t5, first find its position at t4.
v_x1 = v * r_x
v_y1 = v * r_y
v_z1 = v * r_z
x_bird_4 = 4 * v_x1 + 10 * v
y_bird_4 = 4 * v_y1 + v
z_bird_4 = 4 * v_z1 + 6 * v

# Then, calculate the bird's position at t5 by propagating its motion from t4.
v_x5_orig = -v * math.cos(gamma_rad)
v_z5_orig = -v * math.sin(gamma_rad)
x_bird_5 = x_bird_4 + v_x5_orig * (t5 - t4)
y_bird_5 = y_bird_4
z_bird_5 = z_bird_4 + v_z5_orig * (t5 - t4)

# Step 5: Analyze the final leg of the actual journey from t5 to t6.
delta_t_final = t6 - t5

# The bird travels from its position at t5 to the new meeting point (0, y_meet_new, 0)
# at a constant speed v. We can find its velocity components for this leg.
v_xf = -x_bird_5 / delta_t_final
v_zf = -z_bird_5 / delta_t_final
v_yf_sq = v**2 - v_xf**2 - v_zf**2
v_yf = math.sqrt(v_yf_sq)  # Northward motion implies a positive root.

# This allows us to find the new meeting point's y-coordinate.
y_meet_new = y_bird_5 + v_yf * delta_t_final

# Step 6: Calculate the man's required constant acceleration, a3.
# The man must travel from y_man_5 to y_meet_new in the time interval delta_t_final.
delta_y_man = y_meet_new - y_man_5

# We use the kinematic equation: Δy = v_initial * t + 0.5 * a * t^2
# and solve for the acceleration, a3.
numerator = delta_y_man - u * delta_t_final
denominator = 0.5 * delta_t_final**2
a3 = numerator / denominator

# As requested, output the numbers used in the final equation.
print("The man's final acceleration, a_3, is found using the kinematic equation: Δy = v_initial * t + 0.5 * a_3 * t^2")
print(f"The required displacement for the man is Δy = {y_meet_new:.4f} m - {y_man_5:.4f} m = {delta_y_man:.4f} m.")
print(f"The man's initial velocity is v_initial = {u:.4f} m/s, and the time interval is t = {delta_t_final:.4f} s.")
print("Solving for a_3:")
print(f"a_3 = (Δy - v_initial * t) / (0.5 * t^2)")
print(f"a_3 = ({delta_y_man:.4f} - {u:.4f} * {delta_t_final:.4f}) / (0.5 * {delta_t_final:.4f}**2)")
print(f"a_3 = ({numerator:.4f}) / ({denominator:.4f})")
print(f"The man's constant acceleration during this final leg of his journey is: {a3:.4f} m/s^2")
print(f"<<<{a3:.4f}>>>")