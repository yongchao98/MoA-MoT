# Define the parameters for the transistor with example values.
# C_tg: Top gate capacitance per unit area (in Farads per square meter, F/m^2)
C_tg = 0.015

# V_tg: Top gate voltage (in Volts, V)
V_tg = 1.2

# C_bg: Back gate capacitance per unit area (in Farads per square meter, F/m^2)
C_bg = 0.005

# V_bg: Back gate voltage (in Volts, V)
V_bg = 2.0

# Calculate the displacement field (D) using the derived formula.
# D = C_tg * V_tg + C_bg * V_bg
# The unit of D will be Coulombs per square meter (C/m^2).
D_total = C_tg * V_tg + C_bg * V_bg

# Print the final result, showing the equation with the numbers.
# The format D = C_tg * V_tg + C_bg * V_bg is followed.
print(f"The displacement field D is calculated as:")
print(f"D = C_tg * V_tg + C_bg * V_bg")
print(f"D = {C_tg} * {V_tg} + {C_bg} * {V_bg}")
print(f"D = {C_tg * V_tg} + {C_bg * V_bg}")
print(f"D = {D_total} C/m^2")
