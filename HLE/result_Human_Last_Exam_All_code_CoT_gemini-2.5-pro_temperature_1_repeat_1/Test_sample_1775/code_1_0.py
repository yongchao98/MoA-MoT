# Define the parameters for the calculation.
# All units are in the standard SI system (Meters, Volts, Farads, Coulombs).

# Top gate voltage in Volts (V)
V_tg = 1.2

# Back gate voltage in Volts (V)
V_bg = 0.8

# Top gate capacitance per unit area in Farads per square meter (F/m^2)
C_tg = 1.7E-2

# Back gate capacitance per unit area in Farads per square meter (F/m^2)
C_bg = 0.6E-2

# The transistor channel is grounded (V_channel = 0V).
# The displacement field D is given by the superposition of the fields
# from the back gate and the top gate.
# D = D_bg + D_tg
# Let's define the direction from the back gate to the top gate as positive.
# D_bg = C_bg * (V_bg - V_channel) = C_bg * V_bg
# D_tg = -C_tg * (V_tg - V_channel) = -C_tg * V_tg
# Total D = C_bg * V_bg - C_tg * V_tg

# Calculate the displacement field D
displacement_field = C_bg * V_bg - C_tg * V_tg

# Print the final equation with all the values
print("The displacement field (D) is calculated as:")
print(f"D = C_bg * V_bg - C_tg * V_tg")
print(f"D = {C_bg:.2e} F/m^2 * {V_bg} V - {C_tg:.2e} F/m^2 * {V_tg} V")
print(f"D = {C_bg * V_bg:.4e} C/m^2 - {C_tg * V_tg:.4e} C/m^2")
print(f"D = {displacement_field:.4e} C/m^2")

# The final answer in the format <<<value>>>
# print(f"<<<{displacement_field}>>>")