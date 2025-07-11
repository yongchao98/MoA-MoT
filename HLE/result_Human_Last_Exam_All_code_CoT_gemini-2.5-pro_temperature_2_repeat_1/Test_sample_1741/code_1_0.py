import math

# Define the physical constants and initial conditions.
Q_p = 1.602e-19  # Proton charge in Coulombs
m_p = 1.673e-27  # Proton mass in kg
Ex = 5.0         # Electric field in V/m (x-component)
Bz = 5.0e-4      # Magnetic field in Tesla (z-component)
num_revolutions = 3

# Step 1: Calculate the period of one revolution (T).
# The angular frequency is omega_c = (Q * B) / m.
# The period T = 2 * pi / omega_c = (2 * pi * m) / (Q * B).
period = (2 * math.pi * m_p) / (Q_p * Bz)

# Step 2: Calculate the total time for the given number of revolutions.
total_time = num_revolutions * period

# Step 3: Calculate the drift velocity vector (U_d).
# U_d = (E x B) / B^2.
# E is along a_x, B is along a_z. So E x B is along -a_y.
# U_d = -(Ex / Bz) a_y.
U_d_x = 0.0
U_d_y = -Ex / Bz
U_d_z = 0.0

# Step 4: Calculate the final position.
# After an integer number of revolutions, the net displacement from the circular motion is zero.
# Since the proton starts at the origin, its final position is due entirely to the drift velocity.
# r_final = U_d * total_time.
x_final = U_d_x * total_time
y_final = U_d_y * total_time
z_final = U_d_z * total_time

# Step 5: Print the result in the requested format.
# The final equation for the position vector is r = x*a_x + y*a_y + z*a_z.
print("The final position after three complete revolutions is given by the vector:")
print(f"r = ({x_final:.3f}) a_x + ({y_final:.3f}) a_y + ({z_final:.3f}) a_z m")