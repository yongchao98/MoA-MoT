# Define the parameters with example values in SI units.
# You can change these values to match your specific transistor.

# Ctg: Top gate capacitance per unit area (in Farads per square meter, F/m^2)
# Example: 1.5 uF/cm^2 is equivalent to 1.5e-2 F/m^2
Ctg = 1.5e-2

# Vtg: Top gate voltage (in Volts, V)
Vtg = 2.0

# Cbg: Back gate capacitance per unit area (in F/m^2)
# Example: 0.5 uF/cm^2 is equivalent to 0.5e-2 F/m^2
Cbg = 0.5e-2

# Vbg: Back gate voltage (in Volts, V)
Vbg = 5.0

# The displacement field D is the sum of the charge densities induced by the top and back gates.
# D = (Ctg * Vtg) + (Cbg * Vbg)

# Calculate the total displacement field
displacement_field = (Ctg * Vtg) + (Cbg * Vbg)

# Print the final equation with the numbers substituted, as requested.
# The unit of displacement field is Coulombs per square meter (C/m^2).
print("The displacement field (D) is calculated using the formula: D = C_tg * V_tg + C_bg * V_bg")
print("\nSubstituting the given values:")
print(f"D = {Ctg} * {Vtg} + {Cbg} * {Vbg}")
print(f"\nResulting displacement field:")
print(f"D = {displacement_field} C/m^2")