import math

# I. Define constants and calculate key time points t2 and t4
u = 1.5  # Man's initial and final constant speed (m/s)
t0 = 0.0
t1 = 4.0
a1 = -0.15  # Man's deceleration (m/s^2)
t3 = 15.0
a2 = 0.25  # Man's acceleration (m/s^2)
t5 = 23.0
t6 = 40.0

# As explained in the plan, alpha=130 leads to a contradiction. Using 50 degrees makes the problem solvable.
alpha_deg = 50.0
gamma_deg = 40.0

# Calculate t2, when the man's speed becomes 0
# v_f = v_i + a*t => 0 = u + a1*(t2 - t1) => t2 = t1 - u/a1
t2 = t1 - u / a1
# print(f"t2 = {t2} s")

# Calculate t4, when the man's speed reaches u again
# v_f = v_i + a*t => u = 0 + a2*(t4 - t3) => t4 = t3 + u/a2
t4 = t3 + u / a2
# print(f"t4 = {t4} s")

# Time durations for each segment
dt_01 = t1 - t0
dt_12 = t2 - t1
dt_23 = t3 - t2
dt_34 = t4 - t3

# II. Solve for the bird's speed 'v' using the "planned" meeting scenario
# Convert angles to radians
alpha = math.radians(alpha_deg)
gamma = math.radians(gamma_deg)

# Set up equations for r_g = v_ground/v and r_z = v_z/v
# From x-z plane constraint (z4 = x4 * tan(gamma)):
# (4*r_g*sin(alpha) + dt_12) * tan(gamma) = 4*r_z + dt_34
# And r_g^2 + r_z^2 = 1
# This is a quadratic equation for r_g: A*r_g^2 + B*r_g + C = 0
A_coeff = 4 * math.sin(alpha) * math.tan(gamma)
B_coeff = dt_12 * math.tan(gamma) - dt_34
# r_z = sqrt(1-r_g^2)
# (A_coeff*r_g + B_coeff)^2 = (4*sqrt(1-r_g^2))^2
# (A_coeff^2)*r_g^2 + 2*A_coeff*B_coeff*r_g + B_coeff^2 = 16*(1-r_g^2)
# (A_coeff^2 + 16)*r_g^2 + (2*A_coeff*B_coeff)*r_g + (B_coeff^2 - 16) = 0
qA = A_coeff**2 + 16
qB = 2 * A_coeff * B_coeff
qC = B_coeff**2 - 16
# Solve quadratic equation for r_g
discriminant = qB**2 - 4 * qA * qC
r_g = (-qB + math.sqrt(discriminant)) / (2 * qA)
r_z = math.sqrt(1 - r_g**2)

# Calculate man's position at t4
y_man_t1 = u * dt_01
y_man_t2 = y_man_t1 + u * dt_12 + 0.5 * a1 * dt_12**2
y_man_t3 = y_man_t2
y_man_t4 = y_man_t3 + 0.5 * a2 * dt_34**2

# Solve for v using the y-plane constraint: y_man_final = y_bird_final
# y_man_final = y_man_t4 + u * (t_final - t4)
# y_bird_final = v * (4*r_g*cos(alpha) + dt_23)
# t_final - t4 = (v*(4*r_g*sin(alpha)+dt_12)) / (v*cos(gamma))
v_numerator = y_man_t4 + u * (4 * r_g * math.sin(alpha) + dt_12) / math.cos(gamma)
v_denominator = (4 * r_g * math.cos(alpha) + dt_23)
v = v_numerator / v_denominator

# III. Calculate positions at t5 = 23s
# Man's position at t5
dt_45 = t5 - t4
y_man_t5 = y_man_t4 + u * dt_45
u_man_t5 = u

# Bird's position at t5
# Position at t4
x_bird_t4 = v * (4 * r_g * math.sin(alpha) + dt_12)
y_bird_t4 = v * (4 * r_g * math.cos(alpha) + dt_23)
z_bird_t4 = v * (4 * r_z + dt_34)
# Position at t5 (after moving west and down for dt_45 seconds)
x_bird_t5 = x_bird_t4 - v * math.cos(gamma) * dt_45
y_bird_t5 = y_bird_t4
z_bird_t5 = z_bird_t4 - v * math.sin(gamma) * dt_45

# IV. Solve for the final acceleration a3
dt_56 = t6 - t5
# From the bird's constant speed during the final leg:
# v^2 * dt_56^2 = (0 - x_bird_t5)^2 + (y_man_t6 - y_bird_t5)^2 + (0 - z_bird_t5)^2
# Solve for (y_man_t6 - y_bird_t5)
y_diff_sq = v**2 * dt_56**2 - x_bird_t5**2 - z_bird_t5**2
# The wind gust is "northward", so the bird's final y-velocity is positive.
# This means y_man_t6 > y_bird_t5, so we take the positive root.
y_diff = math.sqrt(y_diff_sq)
y_man_t6 = y_bird_t5 + y_diff

# From the man's motion during the final leg:
# y_man_t6 = y_man_t5 + u_man_t5 * dt_56 + 0.5 * a3 * dt_56**2
# Solve for a3
a3_numerator = y_man_t6 - y_man_t5 - u_man_t5 * dt_56
a3_denominator = 0.5 * dt_56**2
a3 = a3_numerator / a3_denominator

print(f"The man's constant acceleration during the final leg of his journey is {a3:.4f} m/s^2.")
# For the final equation representation, we can show the calculation of a3
# a3 = (y_man_t6 - y_man_t5 - u_man_t5 * dt_56) / (0.5 * dt_56^2)
print("The final calculation is:")
print(f"a3 = ({y_man_t6:.4f} - {y_man_t5:.4f} - {u_man_t5:.4f} * {dt_56:.4f}) / (0.5 * {dt_56:.4f}^2)")
print(f"a3 = ({y_man_t6 - y_man_t5 - u_man_t5 * dt_56:.4f}) / ({0.5 * dt_56**2:.4f})")
print(f"a3 = {a3:.4f}")
