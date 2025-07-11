import math

# Step 1: Define all the given constants
u = 1.5  # m/s, man's initial speed
t0 = 0
t1 = 4.0  # s
alpha = 130.0  # degrees, not directly used due to simplification
a1 = -0.15  # m/s^2, man's deceleration
t3 = 15.0  # s
a2 = 0.25  # m/s^2, man's acceleration
gamma = 40.0  # degrees
t5 = 23.0  # s
t6 = 40.0  # s

# Step 2: Analyze the man's motion until t3
# Position at t1
y_t1 = u * t1
u_t1 = u

# Time to stop (t2)
t_decel = u / (-a1)
t2 = t1 + t_decel

# Position and velocity at t2
u_t2 = u_t1 + a1 * (t2 - t1) # This will be 0
y_t2 = y_t1 + u_t1 * (t2 - t1) + 0.5 * a1 * (t2 - t1)**2

# Position and velocity at t3
y_t3 = y_t2
u_t3 = 0

# Step 3: Determine t4 by assuming the man's speed returns to u
# u(t4) = u_t3 + a2 * (t4 - t3) = u
t4 = t3 + (u - u_t3) / a2

# Step 4: Calculate the man's state (position and velocity) at t4 and t5
y_t4 = y_t3 + u_t3 * (t4 - t3) + 0.5 * a2 * (t4 - t3)**2
u_t4 = u_t3 + a2 * (t4 - t3)

# From t4 to t5, the man was supposed to move at a constant velocity u_t4
y_t5 = y_t4 + u_t4 * (t5 - t4)
u_t5 = u_t4

# Step 5: Analyze the bird's planned path to find its properties
# We solve for k_g = v_ground1 / v_total based on the planned rendezvous geometry.
# The simplified equation is (S^2+1)*k_g^2 - (2*K*S)*k_g + (K^2-1) = 0
S = math.sin(math.radians(gamma))
K = 1.5 - 2.5 * math.tan(math.radians(gamma))
a_quad = S**2 + 1
b_quad = -2 * K * S
c_quad = K**2 - 1
# Solve quadratic equation for k_g
discriminant = b_quad**2 - 4 * a_quad * c_quad
k_g = (-b_quad + math.sqrt(discriminant)) / (2 * a_quad) # take positive root

# Now find the planned rendezvous time, t_final
# This simplifies because sin(130) = cos(40)
t_final_delta = 4 * k_g + 10
t_final = t4 + t_final_delta

# Now find the bird's speed 'v'
# y_bird(t4) must equal y_man(t_final) in the planned rendezvous
y_man_t_final = y_t4 + u_t4 * (t_final - t4)
# y_bird(t4) is also given by v * (1 - 4*k_g*cos(130))
# Note: cos(130) = -cos(50)
y_bird_factor = 1 - 4 * k_g * math.cos(math.radians(alpha))
v = y_man_t_final / y_bird_factor

# Step 6: Determine the actual meeting point y_man(t6)
y_bird_t5 = y_man_t_final # bird's y position doesn't change until t5

# Calculate the bird's new vertical velocity component vy_new
# From v^2 = vx_new^2 + vy_new^2 + vz_new^2 and other relations
# we get (vy_new/v)^2 = 1 - ((t_final-t5)/(t6-t5))^2
vy_new_ratio_sq = 1 - ((t_final - t5) / (t6 - t5))**2
vy_new = v * math.sqrt(vy_new_ratio_sq) # positive root because "northward" shift

# Calculate the final meeting position
y_man_t6 = y_bird_t5 + vy_new * (t6 - t5)

# Step 7: Calculate the final acceleration a3
# Use the kinematic equation: y(t6) = y(t5) + u(t5)*dt + 0.5*a3*dt^2
delta_t_final = t6 - t5
a3 = (2 * (y_man_t6 - y_t5 - u_t5 * delta_t_final)) / (delta_t_final**2)

# Print the final calculation step by step
print("The final leg of the man's journey starts at t5 = {} s and ends at t6 = {} s.".format(t5, t6))
print("The duration of this leg is delta_t = t6 - t5 = {:.1f} s.".format(delta_t_final))
print("At t5, the man's position is y_man(t5) = {:.2f} m.".format(y_t5))
print("At t5, the man's speed is u_man(t5) = {:.2f} m/s.".format(u_t5))
print("The man and bird meet at a final position y_man(t6) = {:.2f} m.".format(y_man_t6))
print("\nTo find the man's constant acceleration (a3), we use the kinematic equation:")
print("a3 = (2 * (y_man(t6) - y_man(t5) - u_man(t5) * delta_t)) / delta_t^2")
print("a3 = (2 * ({:.2f} - {:.2f} - {:.2f} * {:.1f})) / {:.1f}^2".format(
    y_man_t6, y_t5, u_t5, delta_t_final, delta_t_final))
final_numerator = 2 * (y_man_t6 - y_t5 - u_t5 * delta_t_final)
final_denominator = delta_t_final**2
print("a3 = {:.2f} / {:.2f}".format(final_numerator, final_denominator))
print("a3 = {:.3f} m/s^2".format(a3))

print(f"\n<<< {a3:.3f} >>>")