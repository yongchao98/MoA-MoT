import math

# --- Step 0: Define the given constants ---
L = 40.0  # Length of the drawbridge in meters
y = 10.0  # Height of the edge above the ground in meters

# --- Step 1: Establish the relationship and find its derivative ---
# The relationship between height y and angle theta is: y = L * cos(theta)
# Differentiating with respect to time t gives the vertical velocity dy/dt:
# dy/dt = -L * sin(theta) * d(theta)/dt

# --- Step 2: Find sin(theta) at the moment y = 10 m ---
# From y = L * cos(theta), we find cos(theta)
cos_theta = y / L
# From sin^2(theta) + cos^2(theta) = 1, we find sin(theta).
# As the bridge is raised, theta is between 0 and pi/2, so sin(theta) is positive.
sin_theta = math.sqrt(1 - cos_theta**2)

# --- Step 3: Use the given value for d(theta)/dt ---
# The rate of change of theta is given as d(theta)/dt = -(3*pi/10) / cos(pi/12)
dtheta_dt_val = -(3 * math.pi / 10) / math.cos(math.pi / 12)

# --- Step 4: Calculate the vertical velocity dy/dt ---
# Substitute the values into the derivative formula.
# dy/dt = -L * sin(theta) * d(theta)/dt
dy_dt = -L * sin_theta * dtheta_dt_val

# --- Final Output ---
# The problem asks to show each number in the final equation.
# The simplified expression for dy/dt is (3 * pi * sqrt(15)) / cos(pi/12).
print("The final equation for the vertical velocity dy/dt is:")
print("dy/dt = (3 * pi * sqrt(15)) / cos(pi/12)")

# Now, we substitute the numerical values for each part of the expression.
numerator_val = 3 * math.pi * math.sqrt(15)
denominator_val = math.cos(math.pi / 12)

print("\nSubstituting the numerical values into the equation:")
print(f"dy/dt = (3 * {math.pi:.5f} * {math.sqrt(15):.5f}) / {denominator_val:.5f}")
print(f"dy/dt = {numerator_val:.5f} / {denominator_val:.5f}")
print(f"dy/dt = {dy_dt:.5f}")

print(f"\nThe vertical velocity of the moving edge is approximately {dy_dt:.4f} m/min.")