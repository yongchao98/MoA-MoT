import numpy as np

# --- Given constants ---
Qp = 1.602e-19  # Charge of a proton in Coulombs
mp = 1.673e-27  # Mass of a proton in kg
B_vec = np.array([0, 0, 5.0e-4])  # Magnetic field in Tesla (T)
E_vec = np.array([5.0, 0, 0])      # Electric field in Volts/meter (V/m)
U0_vec = np.array([2.5e5, 0, 0])   # Initial velocity in meters/second (m/s)
revolutions = 3

# --- Step 1: Calculate the drift velocity (v_d = (E x B) / B^2) ---
E_cross_B = np.cross(E_vec, B_vec)
B_magnitude_sq = np.dot(B_vec, B_vec)
v_drift = E_cross_B / B_magnitude_sq

# --- Step 2: Calculate the period of one revolution (T = 2*pi/omega) ---
# Cyclotron angular frequency omega = q*B/m
B_magnitude = np.linalg.norm(B_vec)
omega = (Qp * B_magnitude) / mp
# Period T
T = 2 * np.pi / omega

# --- Step 3: Calculate the total time for the given number of revolutions ---
t_final = revolutions * T

# --- Step 4: Calculate the final position (R = v_d * t_final) ---
# The initial position is the origin (0,0,0)
R_final = v_drift * t_final

# --- Print the results ---
print("--- Calculation Steps ---")
print(f"Electric Field (E): {E_vec} V/m")
print(f"Magnetic Field (B): {B_vec} T")
print(f"Drift Velocity (v_d = E x B / B^2): {v_drift} m/s")
print(f"Cyclotron Frequency (omega = qB/m): {omega:.4e} rad/s")
print(f"Period of one revolution (T = 2*pi/omega): {T:.4e} s")
print(f"Total time for {revolutions} revolutions (t_final = {revolutions} * T): {t_final:.4e} s")
print("\n--- Final Position Calculation ---")
print(f"The final position vector R is calculated as v_d * t_final.")
print(f"R = ({v_drift[0]:.4f} ax + {v_drift[1]:.4f} ay + {v_drift[2]:.4f} az) m/s * {t_final:.4f} s")

# To satisfy the "output each number in the final equation" requirement
x_final_calc = f"{v_drift[0]:.1f} * {t_final:.4f} = {R_final[0]:.4f}"
y_final_calc = f"{v_drift[1]:.1f} * {t_final:.4f} = {R_final[1]:.4f}"
z_final_calc = f"{v_drift[2]:.1f} * {t_final:.4f} = {R_final[2]:.4f}"

print(f"x_final = {x_final_calc} m")
print(f"y_final = {y_final_calc} m")
print(f"z_final = {z_final_calc} m")


print("\n--- Final Answer ---")
print(f"The final position of the proton after {revolutions} complete revolutions is:")
print(f"R = {R_final[0]:.4f} ax + {R_final[1]:.4f} ay + {R_final[2]:.4f} az meters")
print(f"R = ({R_final[0]:.4f}, {R_final[1]:.4f}, {R_final[2]:.4f}) meters")

# The final result in the required format
final_position_tuple = (round(R_final[0], 4), round(R_final[1], 4), round(R_final[2], 4))
# <<< (0.0, -3.937, 0.0) >>>