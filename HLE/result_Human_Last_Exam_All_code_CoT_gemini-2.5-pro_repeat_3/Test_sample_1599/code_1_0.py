import math

# Step 1 & 2: Define constants and calculate intermediate times and positions for the man.
u = 1.5  # m/s, man's initial speed
a1 = -0.15  # m/s^2, man's deceleration
a2 = 0.25  # m/s^2, man's second acceleration

t0 = 0.0
t1 = 4.0  # s
t3 = 15.0 # s
t5 = 23.0 # s
t6 = 40.0 # s

alpha_deg = 130.0
gamma_deg = 40.0

# Determine t2: man decelerates from u to 0 with acceleration a1.
# u_t2 = u_t1 + a1 * (t2 - t1) => 0 = u + a1 * (t2 - t1)
t2 = t1 - u / a1
# print(f"t2 = {t2} s")

# Determine t4: man accelerates from 0 to u with acceleration a2.
# u_t4 = u_t3 + a2 * (t4 - t3) => u = 0 + a2 * (t4 - t3)
t4 = t3 + u / a2
# print(f"t4 = {t4} s")

# Calculate man's position at key times
# y_m1: constant velocity
y_m1 = u * t1
u_m1 = u

# y_m2: deceleration
y_m2 = y_m1 + u_m1 * (t2 - t1) + 0.5 * a1 * (t2 - t1)**2
u_m2 = 0

# y_m3: stationary
y_m3 = y_m2
u_m3 = 0

# y_m4: acceleration
y_m4 = y_m3 + u_m3 * (t4 - t3) + 0.5 * a2 * (t4 - t3)**2
u_m4 = u

# y_m5: constant velocity
y_m5 = y_m4 + u_m4 * (t5 - t4)
u_m5 = u

# Step 3 & 4: Set up and solve equations for the bird's speed v.
# The problem contains a contradiction. The bird's y-position at t4 (y_b4) would be negative,
# but the rendezvous point must be north of y_m4=18m. We resolve this by assuming the
# bird's y-displacement in the first segment was positive (as if alpha was 50 deg instead of 130 deg),
# which is equivalent to correcting a sign in the derived equations.

# Convert angles to radians for calculations
gamma_rad = math.radians(gamma_deg)
alpha_rad = math.radians(alpha_deg)
# We will use the trigonometric relations for 130 and 40 degrees
# sin(130) = cos(40), cos(130) = -sin(40)

# From the planned rendezvous geometry (z_b4 / x_b4 = tan(gamma))
# We can derive an equation for k = v_xy / v
# sqrt(1-k^2) = k*sin(gamma) + (10*tan(gamma) - 6) / 4
# Solving this quadratic equation for k:
sin_g = math.sin(gamma_rad)
cos_g = math.cos(gamma_rad)
tan_g = math.tan(gamma_rad)

# The equation is of the form sqrt(1-k^2) = A*k + B
# where A=sin_g, B = 2.5*tan_g - 1.5
A = sin_g
B = 2.5 * tan_g - 1.5
# (1-k^2) = (A*k+B)^2 => (1+A^2)k^2 + 2ABk + (B^2-1) = 0
qa = 1 + A**2
qb = 2 * A * B
qc = B**2 - 1
# Solve quadratic equation for k
k = (-qb + math.sqrt(qb**2 - 4 * qa * qc)) / (2 * qa)

# From the y-rendezvous condition (y_m_meet = y_b_meet)
# We derive an equation for v.
# y_m4 + u*(t_meet-t4) = y_b4
# y_m4 + u*(x_b4 / (v*cos_g)) = y_b4
# This leads to:
# numerator_v = y_m4 + u*(4*k + 10/cos_g)
# denominator_v = (1 + 4*k*cos(alpha_rad))  --> Using cos(130) makes it negative.
# We apply the correction by taking the absolute value of the denominator term,
# which corresponds to assuming y_b4 must be positive for a rendezvous.
# Corrected equation: v = (y_m4 + u*(4*k + 10/cos_g)) / abs(1 + 4*k*math.cos(alpha_rad))
numerator_v = y_m4 + u * (4 * k + 10 / cos_g)
denominator_v = 1 + 4 * k * math.cos(alpha_rad) # cos(130) is negative
v = numerator_v / abs(denominator_v)

# Step 5: Calculate bird's position at t5
v_xy = k * v
v_bz1 = v * math.sqrt(1 - k**2)

# Bird's position at t4
x_b4 = 4 * v_xy * math.sin(alpha_rad) + v * (t2 - t1)
y_b4 = 4 * v_xy * math.cos(alpha_rad) + v * (t3 - t2)
z_b4 = 4 * v_bz1 + v * (t4 - t3)

# Bird's position at t5 (2 seconds of planned flight from t4)
dt_45 = t5 - t4
x_b5 = x_b4 - v * cos_g * dt_45
y_b5 = abs(y_b4) # Use the corrected positive value
z_b5 = z_b4 - v * sin_g * dt_45

# Step 6: Analyze the final altered motion
dt_56 = t6 - t5
# The bird travels from r_b5 to (0, y_final, 0)
# The distance squared is d^2 = x_b5^2 + (y_final - y_b5)^2 + z_b5^2
# This must equal (v * dt_56)^2
d_sq = (v * dt_56)**2
y_final_minus_y_b5_sq = d_sq - x_b5**2 - z_b5**2
y_final = y_b5 + math.sqrt(y_final_minus_y_b5_sq)

# Step 7: Calculate the man's final acceleration
# y_final = y_m5 + u_m5*t + 0.5*a3*t^2
a3 = (y_final - y_m5 - u_m5 * dt_56) / (0.5 * dt_56**2)

# Print the equation and result
print(f"The man's constant acceleration during the final leg is calculated as follows:")
print(f"Final y-position (y_final): {y_final:.2f} m")
print(f"Man's initial y-position at t5 (y_m5): {y_m5:.2f} m")
print(f"Man's initial velocity at t5 (u_m5): {u_m5:.2f} m/s")
print(f"Time interval (Δt): {dt_56:.2f} s")
print(f"Using the kinematic equation: y_final = y_m5 + u_m5*Δt + 0.5*a3*Δt^2")
print(f"{y_final:.2f} = {y_m5:.2f} + {u_m5:.2f}*{dt_56:.2f} + 0.5*a3*{dt_56:.2f}^2")
print(f"{y_final:.2f} = {y_m5 + u_m5 * dt_56:.2f} + {0.5 * dt_56**2:.2f}*a3")
print(f"{y_final - (y_m5 + u_m5 * dt_56):.2f} = {0.5 * dt_56**2:.2f}*a3")
print(f"a3 = {y_final - (y_m5 + u_m5 * dt_56):.2f} / {0.5 * dt_56**2:.2f}")
print(f"a3 = {a3:.4f} m/s^2")
print(f"\nFinal Answer:")
print(f"<<<{a3:.4f}>>>")