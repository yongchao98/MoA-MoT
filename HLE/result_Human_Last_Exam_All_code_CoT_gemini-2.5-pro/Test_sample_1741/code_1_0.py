import math

# --- Given Constants ---
# Magnetic field in T
B = 5.0e-4
# Electric field in V/m
E = 5.0
# Proton charge in Coulombs
Qp = 1.602e-19
# Proton mass in kg
mp = 1.673e-27
# Number of revolutions
n_rev = 3

# --- Calculations ---

# 1. Calculate the angular frequency (cyclotron frequency)
# omega = (q * B) / m
omega = (Qp * B) / mp

# 2. Calculate the period of one revolution
# T = 2 * pi / omega
T = (2 * math.pi) / omega

# 3. Calculate the total time for 3 revolutions
t_final = n_rev * T

# 4. Calculate the drift velocity component in the y-direction
# The drift velocity vector is V_d = (E x B) / B^2 = -(E/B) ay
# So, the y-component of the velocity is vy_drift = -E / B
vy_drift = -E / B

# 5. Calculate the final position
# After 3 full revolutions, the displacement from the circular motion is zero.
# The final position is determined solely by the drift velocity over the total time.
# R = V_d * t_final. The displacement is only in the y-direction.
final_x = 0.0
final_y = vy_drift * t_final
final_z = 0.0

# --- Print the Final Equation and Result ---

print("The final position (R) is found by multiplying the drift velocity (V_d) by the total time for 3 revolutions (t).")
print("R = (V_d) * t")
print(f"R = ({vy_drift:.4e} ay m/s) * ({t_final:.4e} s)")
print("\nCalculating the final position vector components:")
print(f"Rx = {final_x:.4f} m")
print(f"Ry = {final_y:.4f} m")
print(f"Rz = {final_z:.4f} m")
print("\nFinal Position Vector:")
print(f"R = {final_x:.4f} ax + {final_y:.4f} ay + {final_z:.4f} az (m)")