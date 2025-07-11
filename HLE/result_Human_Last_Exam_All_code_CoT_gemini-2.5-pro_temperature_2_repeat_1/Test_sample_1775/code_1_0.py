import sys

# Define example values for the calculation.
# You can change these values to match your specific transistor.

# Capacitance per unit area for the top gate (in Farads per square meter, F/m^2)
Ctg = 1.5e-2  # Corresponds to 1.5 uF/cm^2

# Capacitance per unit area for the back gate (in Farads per square meter, F/m^2)
Cbg = 0.5e-2  # Corresponds to 0.5 uF/cm^2

# Voltage applied to the top gate (in Volts, V)
Vtg = 1.2

# Voltage applied to the back gate (in Volts, V)
Vbg = 5.0

# --- Calculation ---

# The total displacement field (D) through the transistor channel is the sum
# of the displacement fields from the top and back gates.
# Formula: D = Ctg * Vtg + Cbg * Vbg

# Calculate the contribution from each gate
D_tg = Ctg * Vtg
D_bg = Cbg * Vbg

# Calculate the total displacement field
D_total = D_tg + D_bg

# --- Output the result ---

print("This script calculates the total displacement field (D) in a dual-gate FET.")
print("The formula used is: D = (Top Gate Capacitance * Top Gate Voltage) + (Back Gate Capacitance * Back Gate Voltage)")
print("\n--- Input Values ---")
print(f"Top Gate Capacitance per Area (Ctg): {Ctg} F/m^2")
print(f"Top Gate Voltage (Vtg):              {Vtg} V")
print(f"Back Gate Capacitance per Area (Cbg): {Cbg} F/m^2")
print(f"Back Gate Voltage (Vbg):              {Vbg} V")

print("\n--- Calculation Steps ---")
print(f"D = (Ctg * Vtg) + (Cbg * Vbg)")
print(f"D = ({Ctg} * {Vtg}) + ({Cbg} * {Vbg})")
print(f"D = {D_tg} + {D_bg}")
print(f"D = {D_total} C/m^2")

print("\n--- Final Answer ---")
print(f"The total displacement field through the transistor is {D_total:.4g} C/m^2.")
sys.stdout.flush()
