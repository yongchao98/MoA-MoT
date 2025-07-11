import math

# Step 1: Define the given constants from the problem.
L = 40  # Length of the drawbridge in meters
y = 10  # Vertical height of the edge in meters

# Step 2: Determine the values of cos(theta) and sin(theta) at the specified height.
# The relationship between height y, length L, and angle theta is y = L * cos(theta).
# From this, we find cos(theta).
cos_theta = y / L
# Using the identity sin^2(theta) + cos^2(theta) = 1, we find sin(theta).
# sin(theta) is positive as theta is an acute angle in this context.
sin_theta = math.sqrt(1 - cos_theta**2)

# Step 3: Calculate the given rate of change of theta, d(theta)/dt.
# The problem provides the formula for d(theta)/dt. The angle in math.cos must be in radians.
# math.pi / 12 is already in radians.
d_theta_dt = -(3 * math.pi / 10) / math.cos(math.pi / 12)

# Step 4: Calculate the vertical velocity (dy/dt) using the related rates formula.
# The formula is dy/dt = -L * sin(theta) * d(theta)/dt.
dy_dt = -L * sin_theta * d_theta_dt

# Step 5: Print the final answer and the values used in the calculation, as requested.
print("To find the vertical velocity (dy/dt), we use the formula: dy/dt = -L * sin(theta) * d(theta)/dt\n")
print("Here are the values used in the final calculation:")
print(f"L (Length of the bridge) = {L}")
print(f"sin(theta) at y={y}m = {sin_theta:.4f}")
print(f"d(theta)/dt (Given rate of angle change) = {d_theta_dt:.4f}\n")

print(f"The calculation is: dy/dt = -({L}) * ({sin_theta:.4f}) * ({d_theta_dt:.4f})")
print(f"The vertical velocity of the moving edge is: {dy_dt:.4f} m/s")
