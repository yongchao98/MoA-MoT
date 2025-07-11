# Define the parameters for the field effect transistor
# Vtg: Top gate voltage in Volts (V)
Vtg = 1.2
# Ctg: Top gate capacitance per unit area in Farads per square meter (F/m^2)
Ctg = 1.5e-2
# Vbg: Back gate voltage in Volts (V)
Vbg = 2.5
# Cbg: Back gate capacitance per unit area in Farads per square meter (F/m^2)
Cbg = 0.5e-2
# The dielectric constant of the transistor material itself (epsilon_s) is not needed for this calculation.

# Calculate the displacement field (D)
# D = Ctg * Vtg + Cbg * Vbg
# The unit of D is Coulombs per square meter (C/m^2)
D = Ctg * Vtg + Cbg * Vbg

# Print the calculation and the result
print("The formula for the displacement field (D) is: D = Ctg * Vtg + Cbg * Vbg")
print(f"D = {Ctg} F/m^2 * {Vtg} V + {Cbg} F/m^2 * {Vbg} V")
print(f"The displacement field through the transistor is: {D} C/m^2")
