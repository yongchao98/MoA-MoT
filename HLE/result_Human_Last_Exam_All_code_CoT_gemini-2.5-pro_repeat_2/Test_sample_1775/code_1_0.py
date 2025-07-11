# Define the parameters for the field effect transistor
# All units are in SI (Meters, Volts, Farads, Coulombs)

# Capacitance per unit area for the top gate (in F/m^2)
# Example: 1.5 uF/cm^2 = 1.5e-6 F / (1e-2 m)^2 = 1.5e-2 F/m^2
C_tg = 1.5e-2

# Capacitance per unit area for the back gate (in F/m^2)
# Example: 0.5 uF/cm^2 = 0.5e-6 F / (1e-2 m)^2 = 0.5e-2 F/m^2
C_bg = 0.5e-2

# Voltage applied to the top gate (in V)
V_tg = 2.0

# Voltage applied to the back gate (in V)
V_bg = -5.0

# The dielectric constant of the transistor (epsilon_s) is not needed for this calculation.

# Calculate the displacement field (D) in the transistor
# D = C_bg * V_bg - C_tg * V_tg
D = C_bg * V_bg - C_tg * V_tg

# Print the final equation and the result
# The unit of displacement field is Coulombs per square meter (C/m^2)
print("The displacement field (D) is calculated by the formula: D = C_bg * V_bg - C_tg * V_tg")
print("Using the given values:")
print(f"D = {C_bg} F/m^2 * {V_bg} V - {C_tg} F/m^2 * {V_tg} V")
print(f"D = {D} C/m^2")
