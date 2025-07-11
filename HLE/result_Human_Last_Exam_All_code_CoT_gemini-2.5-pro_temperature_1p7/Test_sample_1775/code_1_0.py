# Define the physical parameters with example values.
# Ctg: Top gate capacitance per unit area in Farads per square meter (F/m^2)
Ctg = 1.5e-2  # Corresponds to 1.5 uF/cm^2

# Cbg: Back gate capacitance per unit area in Farads per square meter (F/m^2)
Cbg = 0.5e-2  # Corresponds to 0.5 uF/cm^2

# Vtg: Top gate voltage in Volts (V)
Vtg = 1.0

# Vbg: Back gate voltage in Volts (V)
Vbg = -2.0

# The displacement field (D) inside the transistor is the sum of the
# charge densities induced by the top and back gates.
# The formula is D = Ctg * Vtg + Cbg * Vbg.
# The unit of D will be Coulombs per square meter (C/m^2).

# Calculate the displacement field
D = Ctg * Vtg + Cbg * Vbg

# Print the explanation and the final equation with values
print("The displacement field (D) is calculated using the formula: D = Ctg * Vtg + Cbg * Vbg")
print("Using the provided values:")
print(f"D = {Ctg:.2e} F/m^2 * {Vtg:.2f} V + {Cbg:.2e} F/m^2 * {Vbg:.2f} V")

# Print the result
print("\nThe calculated displacement field is:")
print(f"D = {D:.2e} C/m^2")
