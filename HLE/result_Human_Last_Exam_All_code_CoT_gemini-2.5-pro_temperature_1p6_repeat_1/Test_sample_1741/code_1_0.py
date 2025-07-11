import math

# --- Given Constants ---
B_z = 5.0e-4  # Magnetic field strength in Tesla
E_x = 5.0  # Electric field strength in V/m
Q_p = 1.602e-19  # Proton charge in Coulombs
m_p = 1.673e-27  # Proton mass in kg
# Initial velocity U_0 is not needed for final position after integer revolutions

# --- Step 1: Calculate the period of one revolution ---
# The period T depends on the particle's mass, charge, and the magnetic field.
# T = 2 * pi * m / (q * B)
T = (2 * math.pi * m_p) / (Q_p * B_z)

# --- Step 2: Calculate the total time for three revolutions ---
revolutions = 3
t_total = revolutions * T

# --- Step 3: Calculate the drift velocity (U_d) ---
# U_d = (E x B) / B^2
# E = (E_x, 0, 0) and B = (0, 0, B_z)
# E x B = (0 * B_z - 0 * 0, 0 * 0 - E_x * B_z, E_x * 0 - 0 * 0) = (0, -E_x * B_z, 0)
# B^2 = B_z^2
# U_d = (0, -E_x / B_z, 0)
U_d_y = -E_x / B_z
U_d_x = 0
U_d_z = 0

# --- Step 4: Calculate the final position ---
# After 3 complete revolutions, the net displacement is from the drift velocity.
# r_final = U_d * t_total
# Since the particle starts at the origin, the final position is this displacement.
x_final = U_d_x * t_total
y_final = U_d_y * t_total
z_final = U_d_z * t_total

# --- Step 5: Print the results ---
print("--- Calculation Steps ---")
print(f"Period of one revolution (T): {T:.4g} s")
print(f"Total time for {revolutions} revolutions (t_total): {t_total:.4g} s")
print(f"Drift velocity vector (U_d): ({U_d_x:.1f}, {U_d_y:.1f}, {U_d_z:.1f}) m/s")
print("\n--- Final Position Calculation ---")
print("The final position vector r is calculated as r = U_d * t_total.")
print(f"r = ({U_d_x:.1f} ax + {U_d_y:.1f} ay + {U_d_z:.1f} az) m/s * {t_total:.6g} s")

# Final Equation Output
print(f"x = {U_d_x:.1f} m/s * {t_total:.6g} s = {x_final:.3f} m")
print(f"y = {U_d_y:.1f} m/s * {t_total:.6g} s = {y_final:.3f} m")
print(f"z = {U_d_z:.1f} m/s * {t_total:.6g} s = {z_final:.3f} m")

print("\n--- Final Answer ---")
print(f"The final position of the proton after three revolutions is ({x_final:.3f}, {y_final:.3f}, {z_final:.3f}) m.")