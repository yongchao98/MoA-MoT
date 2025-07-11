import numpy as np

# Step 1: Define the given constants and parameters
# Proton properties
Qp = 1.602e-19  # Charge of a proton in Coulombs
mp = 1.673e-27  # Mass of a proton in kg

# Field values
E_vec = np.array([5.0, 0, 0])  # Electric field in V/m (5.0 a_x)
B_vec = np.array([0, 0, 5.0e-4])  # Magnetic field in Tesla (5.0e-4 a_z)
B_mag = np.linalg.norm(B_vec)

# Initial conditions
U0_vec = np.array([2.5e5, 0, 0])  # Initial velocity in m/s (2.5e5 a_x)
r0_vec = np.array([0, 0, 0])      # Initial position (origin)

# Step 2: Calculate the drift velocity
# Ud = (E x B) / B^2
Ud_vec = np.cross(E_vec, B_vec) / (B_mag**2)

# Step 3: Calculate the period of one revolution
# T = 2 * pi * m / (q * B)
T = (2 * np.pi * mp) / (Qp * B_mag)

# Step 4: Calculate the total time for three revolutions
num_revolutions = 3
t_total = num_revolutions * T

# Step 5: Calculate the final position
# The displacement from the gyration after 3 full revolutions is zero.
# The final position is the initial position plus the displacement due to drift.
r_final = r0_vec + Ud_vec * t_total

# Output the results
print("This script calculates the position of a proton after three revolutions in E and B fields.")
print("-" * 50)
print(f"Given values:")
print(f"Proton Charge (Qp): {Qp:.4g} C")
print(f"Proton Mass (mp): {mp:.4g} kg")
print(f"Electric Field (E): {E_vec[0]:.2f} a_x V/m")
print(f"Magnetic Field (B): {B_vec[2]:.2e} a_z T")
print("-" * 50)
print("Calculations:")
print(f"The period of one revolution (T) is calculated as 2*pi*m / (q*B).")
print(f"T = (2 * pi * {mp:.4g}) / ({Qp:.4g} * {B_mag:.2e}) = {T:.4g} s")
print(f"The total time for {num_revolutions} revolutions is t = {num_revolutions} * T = {t_total:.4g} s")
print(f"The drift velocity (Ud) is calculated as (E x B) / B^2.")
print(f"Ud = ({E_vec[0]:.2f} a_x V/m) x ({B_vec[2]:.2e} a_z T) / ({B_mag:.2e})^2 = {Ud_vec[1]:.4g} a_y m/s")
print("-" * 50)
print("Final Position Calculation:")
print("The final position r = r0 + Ud * t.")
print("Since r0 is the origin and the displacement from gyration over 3 revolutions is zero, the position is determined by the drift.")
print(f"r_final_x = r0_x + Ud_x * t = {r0_vec[0]:.3f} + ({Ud_vec[0]:.3f} * {t_total:.4g}) = {r_final[0]:.3f} m")
print(f"r_final_y = r0_y + Ud_y * t = {r0_vec[1]:.3f} + ({Ud_vec[1]:.1f} * {t_total:.4g}) = {r_final[1]:.3f} m")
print(f"r_final_z = r0_z + Ud_z * t = {r0_vec[2]:.3f} + ({Ud_vec[2]:.3f} * {t_total:.4g}) = {r_final[2]:.3f} m")
print("-" * 50)
print(f"The final position vector is: r = {r_final[0]:.3f} a_x + {r_final[1]:.3f} a_y + {r_final[2]:.3f} a_z meters.")