import math

# Define the given constants and initial conditions
Qp = 1.602e-19  # Charge of a proton in Coulombs
mp = 1.673e-27  # Mass of a proton in kg
Bz = 5.0e-4     # Magnetic field in the z-direction in Tesla
Ex = 5.0        # Electric field in the x-direction in V/m
# The initial velocity is not needed to find the position after full revolutions,
# as the net displacement is only due to the E x B drift.

# Step 1: Calculate the period of one revolution (T)
# T = 2 * pi * mp / (Qp * Bz)
T = (2 * math.pi * mp) / (Qp * Bz)

# Step 2: Calculate the total time for three revolutions
num_revolutions = 3
t_final = num_revolutions * T

# Step 3: Calculate the E x B drift velocity (U_d)
# E is in x-dir, B is in z-dir, so E x B is in y-dir.
# The formula is U_d = (E x B) / B^2.
# The direction of E x B is (Ex * ax) x (Bz * az) = Ex * Bz * (ax x az) = -Ex * Bz * ay
# So, U_d = (-Ex * Bz / Bz^2) * ay = -(Ex / Bz) * ay
U_dy = -Ex / Bz
# U_dx and U_dz are zero.

# Step 4: Calculate the final position
# r = U_d * t_final
x_final = 0.0
y_final = U_dy * t_final
z_final = 0.0

# --- Output the result with the calculation breakdown ---
print("The final position of the proton is determined by the E-cross-B drift over the time of three revolutions.")
print("\nThe final position has components (x, y, z). The x and z components are 0.")
print("\nThe final y-position is calculated as:")
print(f"y_final = (Drift Velocity in y) * (Total Time)")
print(f"y_final = -(E_x / B_z) * (3 * 2 * pi * m_p / (Q_p * B_z))")
print(f"y_final = -(6 * pi * m_p * E_x) / (Q_p * B_z^2)\n")

# Print the equation with numerical values substituted
print("Substituting the values:")
print(f"y_final = -(6 * {math.pi:.5f} * {mp:.4e} * {Ex}) / ({Qp:.4e} * ({Bz:.1e})^2)")
print(f"y_final = -({6 * math.pi * mp * Ex}) / ({Qp * Bz**2})")
print(f"y_final = {y_final:.4f} meters\n")

# Print the final position vector
print(f"The final position vector r is ({x_final:.1f}, {y_final:.4f}, {z_final:.1f}) meters.")
print("<<<({:.1f}, {:.4f}, {:.1f})>>>".format(x_final, y_final, z_final))