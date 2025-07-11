import math

# --- Given Information ---
L = 40.0  # Length of the drawbridge in meters
y_target = 10.0  # Target vertical height in meters

# The rate of change of theta is given as a formula.
# d(theta)/dt = - (3 * pi / 10) / cos(pi/12)
# We need the value of its components to calculate it.
cos_pi_div_12 = math.cos(math.pi / 12)
pi_val = math.pi
d_theta_dt = -(3 * pi_val / 10) / cos_pi_div_12

# --- Step 1: Find sin(theta) at the specified height ---
# The vertical height 'y' is related to the angle 'theta' by the equation: y = L * cos(theta)
# We can find cos(theta) when y = 10 m.
cos_theta = y_target / L
# Using the identity sin^2(theta) + cos^2(theta) = 1, we find sin(theta).
# Since the bridge is being raised from horizontal, theta is between pi/2 and 0, so sin(theta) is positive.
sin_theta = math.sqrt(1 - cos_theta**2)

# --- Step 2: Find the vertical velocity dy/dt ---
# The vertical velocity is the time derivative of y.
# dy/dt = d/dt (L * cos(theta))
# Using the chain rule, we get: dy/dt = -L * sin(theta) * d(theta)/dt
# Now we substitute the known values into this equation.
dy_dt = -L * sin_theta * d_theta_dt

# --- Step 3: Print the calculation step-by-step ---
print("The equation for the vertical velocity is: dy/dt = -L * sin(theta) * d(theta)/dt")
print("\nFirst, we find the values for each component at the moment when the bridge is 10 m high:")
print(f"L = {L}")

# Print sin(theta) calculation
print(f"When y = {y_target}, cos(theta) = y/L = {y_target}/{L} = {cos_theta}")
print(f"sin(theta) = sqrt(1 - cos(theta)^2) = sqrt(1 - {cos_theta:.2f}^2) = {sin_theta:.6f}")

# Print d(theta)/dt calculation
print(f"d(theta)/dt is given as - (3 * pi / 10) / cos(pi/12)")
print(f"The value of d(theta)/dt is calculated as: {d_theta_dt:.6f} rad/min")

# Print the final equation with all numbers plugged in
print("\nNow, we substitute these numbers into the velocity equation:")
print("dy/dt = - (L) * (sin(theta)) * (d(theta)/dt)")
print(f"dy/dt = - ({L}) * ({sin_theta:.6f}) * ({d_theta_dt:.6f})")

# Calculate and print the final result
print(f"\nThe calculated vertical velocity is:")
print(f"dy/dt = {dy_dt:.4f} m/min")
