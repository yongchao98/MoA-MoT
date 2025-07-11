# Define the parameters for the field-effect transistor
# V_tg: Top gate voltage in Volts (V)
# V_bg: Bottom gate voltage in Volts (V)
# C_tg: Top gate capacitance per unit area in Farads per square meter (F/m^2)
# C_bg: Bottom gate capacitance per unit area in Farads per square meter (F/m^2)

V_tg = 1.5  # Volts
V_bg = 2.0  # Volts
C_tg = 1.5e-2 # F/m^2 (equivalent to 1.5 uF/cm^2)
C_bg = 1.0e-2 # F/m^2 (equivalent to 1.0 uF/cm^2)

# Calculate the displacement field (D)
# D = C_tg * V_tg + C_bg * V_bg
# The unit of D will be Coulombs per square meter (C/m^2)
displacement_field = C_tg * V_tg + C_bg * V_bg

# Print the final equation with all the numbers
# The numbers are extracted from the variables defined above.
term1 = C_tg * V_tg
term2 = C_bg * V_bg

print(f"The displacement field (D) is calculated as:")
print(f"D = (C_tg * V_tg) + (C_bg * V_bg)")
print(f"D = ({C_tg:.2e} F/m^2 * {V_tg} V) + ({C_bg:.2e} F/m^2 * {V_bg} V)")
print(f"D = {term1} C/m^2 + {term2} C/m^2")
print(f"D = {displacement_field} C/m^2")
