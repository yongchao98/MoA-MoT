import math

# --- Given constants ---
u_initial = 1.5  # m/s
alpha_deg = 50.0 # degrees, interpreted from 130 deg (see plan)
gamma_deg = 40.0 # degrees
a1 = -0.15 # m/s^2
a2 = 0.25 # m/s^2

t0 = 0.0
t1 = 4.0
t3 = 15.0
t5 = 23.0
t6 = 40.0

# --- Plan Step 1: Analyze Man's Journey to find state at t5 ---

# Man's position at t1
y_m1 = u_initial * (t1 - t0)

# Find t2 (when man stops)
# u(t2) = u_initial + a1 * (t2 - t1) = 0
t2 = t1 - u_initial / a1

# Man's position at t2
y_m2 = y_m1 + u_initial * (t2 - t1) + 0.5 * a1 * (t2 - t1)**2
y_m3 = y_m2 # Stays still until t3

# Find t4 (when man reaches u_initial again)
# u(t4) = 0 + a2 * (t4 - t3) = u_initial
t4 = t3 + u_initial / a2

# Man's position at t4
y_m4 = y_m3 + 0.5 * a2 * (t4 - t3)**2

# Man's position and velocity at t5
# From t4 to t5, he moves at constant speed u_initial
dt_4_5 = t5 - t4
y_m5 = y_m4 + u_initial * dt_4_5
u_m5 = u_initial

# --- Plan Step 2 & 3: Analyze Bird's Ideal Journey to find v ---

# Using cos(50) = sin(40) and sin(50) = cos(40) for simplification
C = math.cos(math.radians(gamma_deg)) # cos(40)
S = math.sin(math.radians(gamma_deg)) # sin(40)
T = math.tan(math.radians(gamma_deg)) # tan(40)

# The system of equations for R (=v_1_xy/v) and dt_ideal leads to a quadratic equation for R:
# (16 + 16*S*S) * R^2 + (8*S*(10*T - 6)) * R + ((10*T-6)^2 - 16) = 0
# A simplified version was derived in the thinking process: 22.611 R^2 + 12.29 R - 10.283 = 0
# Let's solve it
# R = (-b +/- sqrt(b^2 - 4ac)) / 2a
a_quad = 16 * C**2 + 16 * S**2
b_quad = 2 * (4 * S) * (10 * T - 6)
c_quad = (10 * T - 6)**2 - 16
# More simply, using derived coefficients:
a_quad_s = 22.611
b_quad_s = 12.29
c_quad_s = -10.283

discriminant = b_quad_s**2 - 4 * a_quad_s * c_quad_s
R = (-b_quad_s + math.sqrt(discriminant)) / (2 * a_quad_s) # Take positive root for speed ratio

# Now find dt_ideal (duration of ideal final leg) and bird speed v
# From dt_ideal * C = 4 * R * C + 10 => dt_ideal = 4*R + 10/C
dt_ideal = 4 * R + 10 / C

# From (18 + 1.5*dt_ideal)/v = 4*R*S + 1
v = (18 + 1.5 * dt_ideal) / (4 * R * S + 1)

# --- Plan Step 4: Determine Actual Meeting Point ---

# Man's position at t4 determines bird's y-position at t4
y_b4 = y_m4 + u_initial * dt_ideal

# Calculate the final meeting y-coordinate, y_m(t6)
# y_m(t6) = y_b4 + v * sqrt( (t6-t5)^2 - (dt_ideal - (t5-t4))^2 )
dt_5_6 = t6 - t5
dt_4_5_bird = t5 - t4
y_m6 = y_b4 + v * math.sqrt(dt_5_6**2 - (dt_ideal - dt_4_5_bird)**2)

# --- Plan Step 5: Calculate Final Acceleration ---

# Use kinematic equation for the man's final leg
# y_m6 = y_m5 + u_m5 * dt_5_6 + 0.5 * a3 * dt_5_6**2
# Solve for a3
a3 = (y_m6 - y_m5 - u_m5 * dt_5_6) / (0.5 * dt_5_6**2)

# --- Output the results ---

print("Step 1: Man's State at t5")
print(f"The man's journey is analyzed to find his state at t5 = {t5} s.")
print(f" - Position at t5: y_m(t5) = {y_m5:.2f} m")
print(f" - Velocity at t5: u_m(t5) = {u_m5:.2f} m/s\n")

print("Step 2: Final Meeting Point")
print(f"The bird's journey is analyzed to find the final meeting point at t6 = {t6} s.")
print(f" - Bird's speed: v = {v:.2f} m/s")
print(f" - Final meeting y-coordinate: y_m(t6) = {y_m6:.2f} m\n")

print("Step 3: Final Acceleration Calculation")
print("The man's final acceleration, a3, is found using the kinematic equation:")
print(f"y_m(t6) = y_m(t5) + u_m(t5) * (t6 - t5) + 0.5 * a3 * (t6 - t5)^2\n")

print("Plugging in the calculated numbers:")
print(f"{y_m6:.2f} = {y_m5:.2f} + {u_m5:.2f} * {dt_5_6} + 0.5 * a3 * {dt_5_6}^2")
print(f"{y_m6:.2f} = {y_m5:.2f} + {u_m5 * dt_5_6:.2f} + {0.5 * dt_5_6**2:.2f} * a3\n")

print(f"Solving for a3:")
print(f"a3 = ({y_m6:.2f} - {y_m5:.2f} - {u_m5 * dt_5_6:.2f}) / {0.5 * dt_5_6**2:.2f}")
print(f"a3 = {y_m6 - y_m5 - u_m5 * dt_5_6:.2f} / {0.5 * dt_5_6**2:.2f}")
print(f"a3 = {a3:.4f} m/s^2")

# Final answer in specified format
print(f"\n<<<{a3:.4f}>>>")