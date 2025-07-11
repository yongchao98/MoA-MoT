import math

# Plan:
# 1. Define the given physical constants for the proton and the fields.
# 2. Calculate the period of one revolution (T) using the formula T = 2*pi*m / (Q*B).
# 3. Calculate the total time for three revolutions (t_final = 3 * T).
# 4. Calculate the drift velocity (U_d) using the formula U_d = -E/B in the y-direction.
# 5. Calculate the final position (r = t_final * U_d), which will only have a y-component.
# 6. Print the final coordinates (x, y, z).

# Step 1: Define constants
B = 5.0e-4      # Magnetic field magnitude in Tesla
E = 5.0         # Electric field magnitude in V/m
Q_p = 1.602e-19 # Proton charge in Coulombs
m_p = 1.673e-27 # Proton mass in kg
revolutions = 3

# Step 2: Calculate the period of one revolution
# The angular frequency is omega = (Q_p * B) / m_p
# The period T = 2 * pi / omega
T = (2 * math.pi * m_p) / (Q_p * B)

# Step 3: Calculate the total time for three revolutions
t_final = revolutions * T

# Step 4: Calculate the drift velocity.
# The drift is in the -y direction with magnitude E/B.
U_d_y = -E / B

# Step 5: Calculate the final position.
# The displacement is due to the drift velocity over the total time.
# The particle starts at the origin.
x_final = 0.0
y_final = t_final * U_d_y
z_final = 0.0

# Step 6: Print the final position coordinates.
# The final equation for the position is r = (x_final, y_final, z_final).
# We output each number (coordinate) of this final position vector.
print("The final position of the proton is given by the coordinates:")
print(f"x = {x_final:.4f} m")
print(f"y = {y_final:.4f} m")
print(f"z = {z_final:.4f} m")