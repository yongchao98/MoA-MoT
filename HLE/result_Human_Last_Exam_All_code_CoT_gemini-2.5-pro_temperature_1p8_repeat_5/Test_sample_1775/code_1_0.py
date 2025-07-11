# Define the parameters for the field-effect transistor.
# These are example values. You can change them to fit a specific scenario.

# C_tg: Top gate capacitance per unit area in Farads per square meter (F/m^2)
C_tg = 1.5e-2  # Example: 1.5 uF/cm^2 = 1.5e-2 F/m^2

# C_bg: Back gate capacitance per unit area in Farads per square meter (F/m^2)
C_bg = 0.5e-2  # Example: 0.5 uF/cm^2 = 0.5e-2 F/m^2

# V_tg: Voltage applied to the top gate in Volts (V)
V_tg = 1.0

# V_bg: Voltage applied to the back gate in Volts (V)
V_bg = 2.0

# The displacement field (D) is given by the formula: D = C_tg * V_tg + C_bg * V_bg
# This is because the total induced charge density on the grounded channel
# is the sum of the charge induced by the top and back gates.

# Calculate the displacement field D
# The unit of D will be Coulombs per square meter (C/m^2)
D = C_tg * V_tg + C_bg * V_bg

# Print the final equation with the numbers substituted in
print("Displacement Field Calculation:")
print(f"D = C_tg * V_tg + C_bg * V_bg")
print(f"D = ({C_tg} F/m^2) * ({V_tg} V) + ({C_bg} F/m^2) * ({V_bg} V)")

# Print the final result
# The unit F*V is a Coulomb (C), so the final unit is C/m^2
print("\nFinal Result:")
print(f"The displacement field (D) is {D} C/m^2")
