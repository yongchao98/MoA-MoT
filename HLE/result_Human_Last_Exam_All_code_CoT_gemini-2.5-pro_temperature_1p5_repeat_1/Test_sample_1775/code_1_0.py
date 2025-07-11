# Description:
# This script calculates the total displacement field (D) in a dual-gate transistor.
# The displacement field is equivalent to the total charge density induced in the channel
# by the top and back gates.

# --- User-definable parameters ---

# Top gate capacitance per unit area (in Farads per square meter, F/m^2)
Ctg = 1.7e-2  # Example: 1.7e-2 F/m^2 for ~2nm SiO2

# Back gate capacitance per unit area (in Farads per square meter, F/m^2)
Cbg = 3.4e-4  # Example: 3.4e-4 F/m^2 for ~100nm SiO2

# Top gate voltage (in Volts, V)
Vtg = 1.5

# Back gate voltage (in Volts, V)
Vbg = 10.0

# --- Calculation ---

# The total displacement field D is the sum of the contributions from each gate.
# Formula: D = Ctg * Vtg + Cbg * Vbg
D = Ctg * Vtg + Cbg * Vbg

# --- Output the result ---

print("Calculation of the Displacement Field (D)")
print("-" * 40)
print(f"The formula is: D = Ctg * Vtg + Cbg * Vbg")
print(f"Plugging in the values:")
# The final equation with each number is printed as requested.
print(f"D = {Ctg} * {Vtg} + {Cbg} * {Vbg}")
print("-" * 40)
print(f"The calculated displacement field is: {D:.4f} C/m^2")