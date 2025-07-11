# Define the given parameters with example values
# Ctg: Top gate capacitance per area in Farads per square meter (F/m^2)
Ctg = 1.7E-2

# Vtg: Top gate voltage in Volts (V)
Vtg = 1.2

# Cbg: Back gate capacitance per area in Farads per square meter (F/m^2)
Cbg = 4.0E-3

# Vbg: Back gate voltage in Volts (V)
Vbg = 10.0

# Calculate the displacement field contributed by each gate
# The unit of displacement field (D) is Coulombs per square meter (C/m^2)
D_tg = Ctg * Vtg
D_bg = Cbg * Vbg

# The total displacement field is the sum of the contributions from both gates
D_total = D_tg + D_bg

# Output the explanation and the final equation with all numbers
print("The total displacement field (D) is the sum of the fields from the top and back gates.")
print("Formula: D = (Ctg * Vtg) + (Cbg * Vbg)")
print("\n--- Calculation with example values ---")
print(f"Top gate contribution: D_tg = {Ctg} F/m^2 * {Vtg} V = {D_tg:.4f} C/m^2")
print(f"Back gate contribution: D_bg = {Cbg} F/m^2 * {Vbg} V = {D_bg:.4f} C/m^2")
print("\nFinal Equation:")
print(f"D = {D_tg:.4f} C/m^2 + {D_bg:.4f} C/m^2")
print(f"D_total = {D_total:.4f} C/m^2")
