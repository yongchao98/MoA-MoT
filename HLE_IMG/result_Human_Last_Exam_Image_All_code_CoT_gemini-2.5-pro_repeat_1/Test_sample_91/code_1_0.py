import math

# Define the given constants and variables
L = 40.0  # Length of the drawbridge in meters
y = 10.0  # Height of the bridge's edge in meters

# Step 1: Establish the geometric relationship
# The vertical height 'y' is given by y = L * cos(θ)

# Step 2: Find the expression for vertical velocity dy/dt
# Differentiating with respect to time t: dy/dt = -L * sin(θ) * dθ/dt

# Step 3: Determine sin(θ) at the given height y = 10 m
# From y = L * cos(θ), we get cos(θ) = y / L
cos_theta = y / L
# Using the identity sin²(θ) + cos²(θ) = 1
# sin(θ) = sqrt(1 - cos²(θ)). We take the positive root as θ is between 0 and π/2.
sin_theta = math.sqrt(1 - cos_theta**2)

# Step 4: Use the given value for dθ/dt
# The problem states dθ/dt = - (3 * π) / (10 * cos(π/12))
dtheta_dt = - (3 * math.pi) / (10 * math.cos(math.pi / 12))

# Step 5: Calculate the vertical velocity dy/dt
dy_dt = -L * sin_theta * dtheta_dt

# Print the calculation steps and the final result
print("The relationship between the vertical height y and the angle θ is: y = L * cos(θ)")
print("The vertical velocity is the time derivative: dy/dt = -L * sin(θ) * dθ/dt\n")

print(f"At the moment when the height y = {y} m:")
print(f"cos(θ) = y / L = {y} / {L} = {cos_theta}")
print(f"sin(θ) = sqrt(1 - {cos_theta}^2) = {sin_theta}\n")

print(f"The given rate of change of the angle is dθ/dt = {dtheta_dt}\n")

print("Substituting these values into the velocity equation:")
print(f"dy/dt = -L * sin(θ) * dθ/dt")
print(f"dy/dt = -({L}) * ({sin_theta}) * ({dtheta_dt})")
print(f"dy/dt = {dy_dt}\n")

print(f"The vertical velocity of the moving edge is {dy_dt} m/s.")
<<<37.79254682710332>>>