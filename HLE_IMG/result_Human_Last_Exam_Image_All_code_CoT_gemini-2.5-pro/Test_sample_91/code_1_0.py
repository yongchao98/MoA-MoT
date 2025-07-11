import math

# Given values
L = 40.0  # Length of the drawbridge in meters
h = 10.0  # Height of the bridge edge above the ground in meters

# 1. Find cos(theta) and sin(theta) at the given height
# h = L * cos(theta) => cos(theta) = h / L
cos_theta = h / L
# sin^2(theta) + cos^2(theta) = 1 => sin(theta) = sqrt(1 - cos^2(theta))
sin_theta = math.sqrt(1 - cos_theta**2)

# 2. Use the given rate of change of theta
# d(theta)/dt = - (3*pi/10) / cos(pi/12)
# We need to calculate the components of this rate
pi = math.pi
cos_pi_12 = math.cos(pi / 12)
d_theta_dt_val = -(3 * pi / 10) / cos_pi_12

# 3. Calculate the vertical velocity
# vh = -L * sin(theta) * d(theta)/dt
vh = -L * sin_theta * d_theta_dt_val

# 4. Print the equation with the numerical values
print("The equation for the vertical velocity (vh) is:")
print("vh = -L * sin(theta) * d(theta)/dt")
print("\nSubstituting the values:")
# Note: sin(theta) = sqrt(1 - (1/4)^2) = sqrt(15)/4
# Note: d(theta)/dt = -(3*pi/10) / cos(pi/12)
print(f"vh = -({L}) * (sqrt(15)/4) * (- (3 * pi / 10) / cos(pi/12))")
print(f"vh = -({L}) * ({sin_theta:.4f}) * ({d_theta_dt_val:.4f})")

# Print the final result
print(f"\nThe vertical velocity of the moving edge of the bridge is {vh:.4f} m/s.")
