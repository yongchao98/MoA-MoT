import math

# --- Given physical constants ---
B_z = 5.0e-4  # Magnetic field in z-direction (T)
E_x = 5.0     # Electric field in x-direction (V/m)
Q_p = 1.602e-19 # Charge of a proton (C)
m_p = 1.673e-27 # Mass of a proton (kg)
num_revolutions = 3

# --- Calculation of the final position ---
# After an integer number of revolutions, the net displacement from the gyration
# is zero. The final position is determined by the E x B drift velocity
# multiplied by the total time for the revolutions.
# Drift velocity in y-direction: U_dy = -E_x / B_z
# Time for n revolutions: t = n * T = n * (2 * pi * m_p) / (Q_p * B_z)
# Final y-position: y = U_dy * t = - (E_x / B_z) * (2 * n * pi * m_p) / (Q_p * B_z)
# This simplifies to y = -(2 * n * pi * m_p * E_x) / (Q_p * B_z^2)

pi = math.pi
x_final = 0.0
z_final = 0.0

# Numerator of the equation for y
numerator_val = -1 * (2 * num_revolutions * pi * m_p * E_x)

# Denominator of the equation for y
denominator_val = Q_p * (B_z**2)

# Final y-position
y_final = numerator_val / denominator_val

# --- Output the results step-by-step ---
print("The final position (x, y, z) is determined by the E x B drift over 3 revolutions.")
print("The x and z components of the final position are 0.\n")
print("The y-component of the final position is given by the equation:")
print("y = - (2 * n * pi * m_p * E_x) / (Q_p * B_z^2)\n")

print("Substituting the numerical values:")
print(f"y = - (2 * {num_revolutions} * {pi:.5f} * {m_p:.4e} kg * {E_x:.1f} V/m) / ({Q_p:.4e} C * ({B_z:.1e} T)^2)\n")

print("Calculating the numerator and denominator separately:")
print(f"Numerator = {numerator_val:.4e}")
print(f"Denominator = {denominator_val:.4e}\n")

print("The final y-position is:")
print(f"y = {numerator_val:.4e} / {denominator_val:.4e} = {y_final:.4f} meters\n")

print("="*40)
print("The final position of the proton is:")
print(f"r = ({x_final:.1f} a_x) + ({y_final:.4f} a_y) + ({z_final:.1f} a_z) meters")
print("="*40)