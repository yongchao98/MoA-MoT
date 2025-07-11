import math

# Step 1: Define constants and initial values
u_initial = 1.5  # m/s
t0 = 0
t1 = 4
a1 = -0.15
t3 = 15
a2 = 0.25
gamma_deg = 40
t5 = 23
t6 = 40

# --- Man's Journey Analysis (t0 to t5) ---

# Segment 1: t0=0 to t1=4 (constant velocity)
y_man_t1 = u_initial * (t1 - t0)
u_man_t1 = u_initial

# Segment 2: t1=4 to t2 (deceleration)
# Calculate t2, the time when the man stops
t2 = t1 - (0 - u_man_t1) / a1
y_man_t2 = y_man_t1 + u_man_t1 * (t2 - t1) + 0.5 * a1 * (t2 - t1)**2
u_man_t2 = 0

# Segment 3: t2 to t3=15 (stationary)
y_man_t3 = y_man_t2
u_man_t3 = 0

# Segment 4: t3=15 to t4 (acceleration)
# Calculate t4, the time when man reaches initial speed again
u_target = u_initial
t4 = t3 + (u_target - u_man_t3) / a2
y_man_t4 = y_man_t3 + u_man_t3 * (t4 - t3) + 0.5 * a2 * (t4 - t3)**2
u_man_t4 = u_target

# Segment 5: t4 to t5=23 (constant velocity)
y_man_t5 = y_man_t4 + u_man_t4 * (t5 - t4)
u_man_t5 = u_man_t4

# --- Bird's Journey Analysis ---

# Step 2: Resolve angle contradiction and define angles in radians
# The problem's description of the bird's first movement is contradictory.
# "upward, eastward, and northward" implies an angle between 0 and 90 degrees with North.
# The given alpha = 130 degrees is used, assuming it is measured from South.
# Angle with North = 180 - 130 = 50 degrees.
alpha_corrected_deg = 50
alpha_rad = math.radians(alpha_corrected_deg)
gamma_rad = math.radians(gamma_deg)

# Step 3: Solve for bird's speed 'v' using the planned rendezvous
# This involves finding k = v_ground_speed / v_airspeed
# by setting up and solving a quadratic equation derived from geometric constraints.
# (4*k*sin(alpha) + 10)*tan(gamma) = 4*sqrt(1-k^2) + (t4-t3)
# Let's solve (4*k*sin(alpha) + 10)*tan(gamma) = 4*sqrt(1-k^2) + 6
# ( (4*sin(alpha)*tan(gamma))k + (10*tan(gamma)-6) )^2 = 16*(1-k^2)
# (2.571k + 2.392)^2 = 16(1-k^2) -> 22.61k^2 + 12.29k - 10.29 = 0
qa = 22.61
qb = 12.29
qc = -10.29
k = (-qb + math.sqrt(qb**2 - 4 * qa * qc)) / (2 * qa) # ground_speed/air_speed ratio

# Calculate time duration of planned final leg
T_planned = (4 * k * math.sin(alpha_rad) + (t2 - t1)) / math.cos(gamma_rad)

# Calculate bird's y-position at t4 in terms of v
y_bird_t4_coeff = (4 * k * math.cos(alpha_rad) + (t3 - t2))

# Calculate bird's y-position at t4 from man's path
y_bird_t4_abs = y_man_t4 + u_man_t4 * T_planned

# Solve for v
v = y_bird_t4_abs / y_bird_t4_coeff

# Step 4: Determine the new meeting point after the gust at t5
# Calculate bird's 3D position at t4
x_bird_t4 = v * (4 * k * math.sin(alpha_rad) + (t2-t1))
y_bird_t4 = y_bird_t4_abs
z_bird_t4 = v * (4 * math.sqrt(1 - k**2) + (t4-t3))

# Calculate bird's 3D position at t5
dt_45 = t5 - t4
x_bird_t5 = x_bird_t4 - v * math.cos(gamma_rad) * dt_45
y_bird_t5 = y_bird_t4 # Bird's northward velocity was zero in this planned segment
z_bird_t5 = z_bird_t4 - v * math.sin(gamma_rad) * dt_45

# Calculate the new meeting point on the y-axis
dt_56 = t6 - t5
# v^2 = v_nx^2 + v_ny^2 + v_nz^2
# (v*dt_56)^2 = x_b5^2 + y_displacement^2 + z_b5^2
y_displacement_from_t5_sq = (v * dt_56)**2 - x_bird_t5**2 - z_bird_t5**2
y_displacement_from_t5 = math.sqrt(y_displacement_from_t5_sq)

# The gust is "northward", so the bird's y-position increases.
y_meet = y_bird_t5 + y_displacement_from_t5

# Step 5: Calculate the man's final acceleration a3
# y_meet = y_man_t5 + u_man_t5 * dt_56 + 0.5 * a3 * dt_56**2
# Rearrange to solve for a3
a3 = (y_meet - y_man_t5 - u_man_t5 * dt_56) * 2 / (dt_56**2)

# Output the final equation with calculated numbers
term1 = y_man_t5
term2 = u_man_t5
term3 = dt_56
term4 = 0.5
term5 = dt_56

print("The final kinematic equation for the man is:")
print(f"y_meet = y_man(t5) + u_man(t5)*(t6-t5) + 0.5*a3*(t6-t5)^2")
print("Plugging in the calculated values:")
print(f"{y_meet:.2f} = {term1:.2f} + {term2:.2f} * {term3:.2f} + {term4:.1f} * a3 * {term5:.2f}^2")
print("\nSolving for the man's final acceleration (a3):")
print(f"a3 = {a3:.3f} m/s^2")