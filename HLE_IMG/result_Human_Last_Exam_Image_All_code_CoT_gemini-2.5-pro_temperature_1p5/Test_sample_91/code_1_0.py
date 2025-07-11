import math

# --- Given values from the problem ---
L = 40  # Length of the drawbridge in meters
h = 10  # Vertical height of the edge in meters

# --- Step 1: Find sin(theta) at the given height ---
# The vertical height h is related to the angle theta by h = L * cos(theta).
# So, cos(theta) = h / L
cos_theta = h / L

# Using the identity sin^2(theta) + cos^2(theta) = 1, we find sin(theta).
# sin(theta) = sqrt(1 - cos^2(theta)). We take the positive root as theta is between 0 and pi/2.
sin_theta = math.sqrt(1 - cos_theta**2)

# --- Step 2: Use the given rate of change of theta ---
# The rate d(theta)/dt is given in the problem statement.
# d(theta)/dt = - (3 * pi / 10) / cos(pi/12)
cos_pi_over_12 = math.cos(math.pi / 12)
dtheta_dt = -(3 * math.pi) / (10 * cos_pi_over_12)

# --- Step 3: Calculate the vertical velocity dh/dt ---
# The vertical velocity dh/dt is the derivative of h = L * cos(theta) with respect to time t.
# Using the chain rule, dh/dt = -L * sin(theta) * d(theta)/dt
vertical_velocity = -L * sin_theta * dtheta_dt

# --- Step 4: Print the explanation and the final answer ---
print("The vertical velocity of the drawbridge's edge (dh/dt) is calculated using the formula derived from h = L * cos(θ):")
print("dh/dt = -L * sin(θ) * dθ/dt")
print("\nFirst, we find the values of the terms for when the height h is 10 m:")
print(f"Given L = {L} m and h = {h} m.")
print(f"cos(θ) = h / L = {h} / {L} = {cos_theta}")
# We calculate sin(theta) from cos(theta)
sin_theta_val = math.sqrt(1 - cos_theta**2)
print(f"sin(θ) = sqrt(1 - {cos_theta}²) = {sin_theta_val:.4f}")

print("\nNext, we calculate the given rate of change for θ:")
print(f"dθ/dt = -3*pi / (10 * cos(pi/12)) = {dtheta_dt:.4f} rad/s")

print("\nFinally, we plug all the numerical values into the velocity equation:")
print(f"dh/dt = -({L}) * ({sin_theta_val:.4f}) * ({dtheta_dt:.4f})")
print(f"dh/dt = {vertical_velocity:.2f} m/s")

# The final answer in the required format
final_answer = round(vertical_velocity, 2)
# The final result is presented in the specified format below.
# <<<37.79>>>