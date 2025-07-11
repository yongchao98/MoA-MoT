import math

# --- Problem Parameters ---
# Length of the drawbridge (m)
L = 40.0
# Target vertical height of the bridge's edge (m)
h_target = 10.0

# --- Derivation and Calculation ---

# Step 1: Define the relationship between height (h), length (L), and angle (theta)
# From the problem's diagram, h = L * cos(theta).

# Step 2: Differentiate with respect to time to find the vertical velocity (v_y = dh/dt).
# dh/dt = d/dt[L * cos(theta)] = -L * sin(theta) * d(theta)/dt

print("The formula for vertical velocity (v_y) is: v_y = -L * sin(theta) * d(theta)/dt")
print("-" * 30)

# Step 3: Find the value of sin(theta) at the moment when h = 10 m.
cos_theta = h_target / L
# sin^2(theta) + cos^2(theta) = 1  => sin(theta) = sqrt(1 - cos^2(theta))
# Angle is between 0 and 90 degrees, so sin(theta) is positive.
sin_theta = math.sqrt(1 - cos_theta**2)

print(f"Given L = {L} m.")
print(f"At height h = {h_target} m, cos(theta) = {h_target}/{L} = {cos_theta}")
print(f"Therefore, sin(theta) = sqrt(1 - {cos_theta}^2) = sqrt(15)/4 ≈ {sin_theta:.4f}")
print("-" * 30)

# Step 4: Use the given value for the rate of change of the angle, d(theta)/dt.
# The expression is given as: -(3*pi/10) / cos(pi/12)
d_theta_dt_numerator = -3 * math.pi
d_theta_dt_denominator = 10 * math.cos(math.pi / 12)
d_theta_dt = (d_theta_dt_numerator / d_theta_dt_denominator)

print("The rate of change of the angle is given as d(theta)/dt = -(3*pi/10) / cos(pi/12)")
print(f"d(theta)/dt ≈ {d_theta_dt:.4f} rad/s")
print("-" * 30)

# Step 5: Substitute all values into the velocity equation.
vertical_velocity = -L * sin_theta * d_theta_dt

print("Substituting the values into the velocity formula:")
# Per instructions, showing the numbers in the final equation.
# Simplified form of the final equation: v_y = (3 * pi * sqrt(15)) / cos(pi/12)
num1 = 3
num2 = 15
num3 = 12
print(f"v_y = -({L}) * (sqrt(1 - {cos_theta}**2)) * (-(3*pi/10) / cos(pi/12))")
print(f"v_y = ({num1} * pi * sqrt({num2})) / cos(pi/{num3})")
print("\nFinal Calculation:")
print(f"v_y = (3 * {math.pi:.4f} * {math.sqrt(15):.4f}) / {math.cos(math.pi/12):.4f}")
print(f"v_y = {vertical_velocity:.4f} m/s")
