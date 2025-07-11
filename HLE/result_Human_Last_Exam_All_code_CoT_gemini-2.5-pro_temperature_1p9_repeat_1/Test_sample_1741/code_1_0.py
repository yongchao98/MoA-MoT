import math

# Define the physical constants and initial conditions
Q_p = 1.602e-19  # Charge of a proton in Coulombs
m_p = 1.673e-27  # Mass of a proton in kg
B_z = 5.0e-4     # Magnetic field in the z-direction in Tesla
E_x = 5.0        # Electric field in the x-direction in V/m
revolutions = 3

# The motion is a superposition of a drift and a circular motion.
# After an integer number of revolutions, the displacement from the circular part is zero.
# The final position is determined by the drift velocity multiplied by the total time.
#
# Final position vector r = (x_final, y_final, z_final)
# x_final = 0
# z_final = 0
# y_final = U_dy * t_total
#
# U_dy (drift velocity in y) = -E_x / B_z
# t_total (time for 3 revolutions) = 3 * T = 3 * (2 * pi / omega_c)
# omega_c (cyclotron frequency) = (Q_p * B_z) / m_p
#
# So, y_final = (-E_x / B_z) * (3 * 2 * pi * m_p) / (Q_p * B_z)
# y_final = - (6 * pi * E_x * m_p) / (Q_p * B_z^2)

# Calculate the final y-position using the derived formula
y_final = - (revolutions * 2 * math.pi * E_x * m_p) / (Q_p * B_z**2)
x_final = 0.0
z_final = 0.0

print("The final position of the proton is calculated based on its drift over three revolution periods.")
print("The final y-coordinate is given by the equation:")
print(f"y_final = - ({revolutions} * 2 * pi * E_x * m_p) / (Q_p * B_z^2)\n")

print("Plugging in the values:")
print(f"y_final = - ({revolutions} * 2 * {math.pi:.5f} * {E_x} * {m_p}) / ({Q_p} * {B_z}**2)")
print(f"y_final = {y_final:.4f} m\n")

print("The final position vector r is:")
print(f"r = ({x_final:.1f}) a_x + ({y_final:.4f}) a_y + ({z_final:.1f}) a_z  (in meters)")