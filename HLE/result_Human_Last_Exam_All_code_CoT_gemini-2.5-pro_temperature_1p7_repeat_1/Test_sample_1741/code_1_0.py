import math

# --- Given physical constants and parameters ---
Qp = 1.602e-19  # Charge of a proton in Coulombs
mp = 1.673e-27  # Mass of a proton in kg
Bz = 5.0e-4     # Magnetic field strength in Tesla (z-direction)
Ex = 5.0        # Electric field strength in V/m (x-direction)
N = 3           # Number of complete revolutions

# --- Final Position Calculation ---
# After N complete revolutions, the displacement from the circular motion is zero.
# The final position is determined by the drift velocity multiplied by the total time.
# x_final and z_final are 0 due to the problem's geometry.
x_final = 0.0
z_final = 0.0

# The final y-position can be calculated using the combined formula:
# y = (Drift Velocity) * (Total Time)
# y = (-Ex / Bz) * (N * Period)
# y = (-Ex / Bz) * (N * 2 * pi / omega_c)
# y = (-Ex / Bz) * (N * 2 * pi * mp / (Qp * Bz))
# This simplifies to the expression used below.

y_final = -(N * 2 * math.pi * mp * Ex) / (Qp * Bz**2)

# --- Print the Result and Equation Steps ---
print("The final position (x, y, z) is determined as follows:")

# Explanation for x and z positions
print(f"x_final = {x_final:.1f} m (Drift is only in y-direction, and gyration displacement is zero after {N} cycles)")
print(f"z_final = {z_final:.1f} m (No forces or initial velocity in the z-direction)")

# Detailed calculation for the y position
print("\ny_final is the total distance traveled at the drift velocity:")
print(f"y_final = - (N * 2 * pi * m_p * E_x) / (Q_p * B_z^2)")
print("Plugging in the numbers:")
print(f"y_final = - ({N} * 2 * {math.pi:.5f} * {mp:.4e} * {Ex:.1f}) / ({Qp:.4e} * ({Bz:.1e})^2)")
print(f"y_final = {y_final:.3f} m")

print(f"\nTherefore, the final position after {N} revolutions is ({x_final:.1f}, {y_final:.3f}, {z_final:.1f}) m.")