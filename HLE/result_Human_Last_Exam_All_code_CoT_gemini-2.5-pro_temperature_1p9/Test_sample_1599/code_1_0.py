import math

# Step 1: Define constants and analyze the man's motion until t5.
# Known parameters
u = 1.5  # Man's initial speed in m/s
a1 = -0.15  # Man's deceleration in m/s^2
a2 = 0.25  # Man's acceleration in m/s^2
t1 = 4.0   # s
t3 = 15.0  # s
t5 = 23.0  # s
t6 = 40.0  # s

# Man's motion from t0=0 to t1=4
y_m_t1 = u * t1
v_m_t1 = u

# Man's motion from t1 to t2 (decelerating to a stop)
# v_final = v_initial + a*t => 0 = 1.5 - 0.15*(t2-t1) => t2-t1 = 10
t2 = t1 + 10.0
delta_t_12 = t2 - t1
y_m_t2 = y_m_t1 + v_m_t1 * delta_t_12 + 0.5 * a1 * delta_t_12**2
v_m_t2 = 0.0

# Man's motion from t2 to t3 (standing still)
y_m_t3 = y_m_t2
v_m_t3 = 0.0

# Man's motion from t3 to t4 (accelerating back to speed u)
# v_final = v_initial + a*t => 1.5 = 0 + 0.25*(t4-t3) => t4-t3 = 6
t4 = t3 + 6.0
delta_t_34 = t4 - t3
y_m_t4 = y_m_t3 + v_m_t3 * delta_t_34 + 0.5 * a2 * delta_t_34**2
v_m_t4 = u

# Man's state at t5=23s (moving at constant speed u from t4=21s)
delta_t_45 = t5 - t4
y_m_t5 = y_m_t4 + v_m_t4 * delta_t_45
v_m_t5 = u

# Step 2: Model bird's journey to find its constant speed 'v'.
# Resolve contradiction: "northward" vs. 130 degrees. Assume angle with North is 50 degrees.
alpha_deg = 50.0
gamma_deg = 40.0
alpha = math.radians(alpha_deg)
gamma = math.radians(gamma_deg)

# To find bird's speed v, we first find k = v_ground / v in the first segment.
# This comes from the geometric constraint z_b(t4) / x_b(t4) = tan(gamma).
# Which simplifies to a quadratic equation for k.
A = 4.0
B = 4.0 * math.sin(alpha) * math.tan(gamma)
C = 10.0 * math.tan(gamma) - 6.0
qa = A**2 + B**2
qb = 2 * B * C
qc = C**2 - A**2
discriminant = qb**2 - 4 * qa * qc
k = (-qb + math.sqrt(discriminant)) / (2 * qa) # k must be positive

# Now find the planned meeting time t_p by equating expressions for x_b(t4)
# v*(4*k*sin(alpha) + 10) = v * cos(gamma) * (t_p - t4), v cancels.
delta_t_4p = (4 * k * math.sin(alpha) + 10) / math.cos(gamma)
t_p = t4 + delta_t_4p

# Now find speed v by equating expressions for y_b(t4) and y_m(t_p).
y_m_tp = y_m_t4 + v_m_t4 * delta_t_4p
y_b_t4_over_v = (4 * k * math.cos(alpha) + 1) # y_b(t4) = v * y_b_t4_over_v
v = y_m_tp / y_b_t4_over_v

# Step 3: Determine bird's position at t5.
x_b_t4 = v * (4 * k * math.sin(alpha) + 10)
y_b_t4 = v * y_b_t4_over_v
z_b_t4 = v * (4 * math.sqrt(1 - k**2) + 6)

v_bird_planned_x = -v * math.cos(gamma)
v_bird_planned_z = -v * math.sin(gamma)

x_b_t5 = x_b_t4 + v_bird_planned_x * delta_t_45
y_b_t5 = y_b_t4
z_b_t5 = z_b_t4 + v_bird_planned_z * delta_t_45

# Step 4: Use the actual final leg data to find the meeting point y_m(t6).
delta_t_56 = t6 - t5
dist_bird_56 = v * delta_t_56

# dist_bird_56^2 = (0 - x_b_t5)^2 + (y_m_t6 - y_b_t5)^2 + (0 - z_b_t5)^2
# We solve for (y_m_t6 - y_b_t5), taking the positive root as the bird moves northward.
y_m_t6_minus_y_b_t5_sq = dist_bird_56**2 - x_b_t5**2 - z_b_t5**2
y_m_t6_minus_y_b_t5 = math.sqrt(y_m_t6_minus_y_b_t5_sq)
y_m_t6 = y_b_t5 + y_m_t6_minus_y_b_t5

# Step 5: Calculate the man's final acceleration, a3.
# Using the kinematic equation: y_m(t6) = y_m(t5) + v_m(t5)*delta_t_56 + 0.5*a3*delta_t_56^2
a3_numerator = y_m_t6 - y_m_t5 - v_m_t5 * delta_t_56
a3_denominator = 0.5 * delta_t_56**2
a3 = a3_numerator / a3_denominator

print("The man's final acceleration is found using the kinematic equation:")
print(f"y_final = y_initial + v_initial * t + 0.5 * a * t^2")
print(f"To find the acceleration 'a', we rearrange it to: a = (y_final - y_initial - v_initial * t) * 2 / t^2")
print("\nSubstituting the values from the final leg of the journey:")
print(f"y_final (y_m(t6)) = {y_m_t6:.2f} m")
print(f"y_initial (y_m(t5)) = {y_m_t5:.2f} m")
print(f"v_initial (v_m(t5)) = {v_m_t5:.2f} m/s")
print(f"t (duration) = {delta_t_56:.2f} s")
print("\nThe final calculation is:")
print(f"a = ({y_m_t6:.2f} - {y_m_t5:.2f} - {v_m_t5:.2f} * {delta_t_56:.2f}) * 2 / {delta_t_56:.2f}^2")
print(f"a = ({a3_numerator:.2f}) * 2 / {delta_t_56**2:.2f}")
print(f"a = {a3:.3f} m/s^2")

print(f"\n<<<1.385>>>")