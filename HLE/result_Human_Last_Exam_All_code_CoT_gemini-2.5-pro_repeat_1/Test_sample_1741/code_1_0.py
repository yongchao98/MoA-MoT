import math

# --- Given Constants ---
# Electric field in V/m
E_x = 5.0
# Magnetic field in Tesla
B_z = 5.0e-4
# Proton charge in Coulombs
Q_p = 1.602e-19
# Proton mass in kg
m_p = 1.673e-27
# Number of revolutions
n = 3

# --- Step 1: Calculate the period of one revolution (T) ---
# The period of gyration in a magnetic field is T = 2 * pi * m / (Q * B)
T = (2 * math.pi * m_p) / (Q_p * B_z)

# --- Step 2: Calculate the total time for n revolutions ---
t_final = n * T

# --- Step 3: Calculate the drift velocity (U_d) ---
# The drift velocity is in the y-direction, given by U_dy = -E_x / B_z
U_dy = -E_x / B_z

# --- Step 4: Calculate the final position ---
# After an integer number of revolutions, the net displacement is from the drift.
# Final position P = U_d * t_final.
# The initial velocity does not affect the final position after complete cycles.
x_final = 0.0
z_final = 0.0
y_final = U_dy * t_final

# --- Step 5: Print the final answer ---
# The final equation for the position is P = x_final*ax + y_final*ay + z_final*az
print("The final position of the proton is given by the equation:")
print(f"Position = ({x_final:.3f} ax) + ({y_final:.3f} ay) + ({z_final:.3f} az) m")
