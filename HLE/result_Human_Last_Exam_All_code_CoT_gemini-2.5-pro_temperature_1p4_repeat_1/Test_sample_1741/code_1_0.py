import math

# Define the given constants
Qp = 1.602e-19  # Proton charge in Coulombs
mp = 1.673e-27  # Proton mass in kg
Bz = 5.0e-4      # Magnetic field in Tesla (in z-direction)
Ex = 5.0         # Electric field in V/m (in x-direction)
U0x = 2.5e5      # Initial velocity in m/s (in x-direction)

print("Step 1: Calculate the period of one revolution (T)")
# The formula for the period is T = 2 * pi * m / (Q * B)
T = (2 * math.pi * mp) / (Qp * Bz)
print(f"T = (2 * pi * {mp:.4g}) / ({Qp:.4g} * {Bz:.1e})")
print(f"T = {T:.4g} s")
print("-" * 30)

print("Step 2: Calculate the total time for three revolutions (t_final)")
num_revolutions = 3
t_final = num_revolutions * T
print(f"t_final = {num_revolutions} * T = {num_revolutions} * {T:.4g}")
print(f"t_final = {t_final:.4g} s")
print("-" * 30)

print("Step 3: Calculate the drift velocity (Ud)")
# The drift velocity vector is Ud = (E x B) / B^2
# E is in x-direction, B is in z-direction, so E x B is in the -y direction.
# Ud_y = -Ex / Bz
Ud_y = -Ex / Bz
print("The drift velocity is in the y-direction, given by Ud_y = -Ex / Bz")
print(f"Ud_y = -{Ex} / {Bz:.1e}")
print(f"Ud_y = {Ud_y:.4g} m/s")
print("-" * 30)

print("Step 4: Calculate the final position (r_final)")
# After an integer number of revolutions, the displacement due to the circular motion is zero.
# The final position is determined solely by the drift.
# r_final = Ud * t_final
final_x = 0.0
final_y = Ud_y * t_final
final_z = 0.0
print("The final position is (x, y, z) where x=0, z=0, and y = Ud_y * t_final")
print(f"y_final = {Ud_y:.4g} * {t_final:.4g}")
print(f"y_final = {final_y:.4g} m")
print("-" * 30)

print("Final Position:")
print(f"The position after three complete revolutions is ({final_x:.1f}, {final_y:.4f}, {final_z:.1f}) m.")
<<<({final_x}, {final_y}, {final_z})>>>