# Define example parameters for the Field Effect Transistor.
# You can change these values to match your specific device.

# V_tg: Top gate voltage in Volts (V)
V_tg = 1.2

# V_bg: Back gate voltage in Volts (V)
V_bg = 0.8

# C_tg: Top gate capacitance per unit area in Farads per square meter (F/m^2)
C_tg = 0.025

# C_bg: Back gate capacitance per unit area in Farads per square meter (F/m^2)
C_bg = 0.008

# --- Calculation ---
# The displacement field (D) at the semiconductor channel surface is the sum
# of the contributions from the top and back gates.
# Formula: D = C_tg * V_tg + C_bg * V_bg
# The unit of D will be Coulombs per square meter (C/m^2), as F*V = C.
# The dielectric constant of the transistor material (epsilon_s) is not needed
# to find D, which is determined by the external gate parameters.

displacement_field = C_tg * V_tg + C_bg * V_bg

# --- Output ---
print("To find the displacement field (D) in the transistor, we use the following formula based on the superposition of the fields from the top and back gates:")
print("D = C_tg * V_tg + C_bg * V_bg")
print("\n--- Using the provided values ---")
print(f"Top Gate Voltage (V_tg): {V_tg} V")
print(f"Back Gate Voltage (V_bg): {V_bg} V")
print(f"Top Gate Capacitance/Area (C_tg): {C_tg} F/m^2")
print(f"Back Gate Capacitance/Area (C_bg): {C_bg} F/m^2")
print("\n--- Step-by-step Calculation ---")
print(f"D = ({C_tg}) * ({V_tg}) + ({C_bg}) * ({V_bg})")
term1 = C_tg * V_tg
term2 = C_bg * V_bg
print(f"D = {term1} + {term2}")
print(f"\nThe final displacement field is:")
print(f"D = {displacement_field} C/m^2")