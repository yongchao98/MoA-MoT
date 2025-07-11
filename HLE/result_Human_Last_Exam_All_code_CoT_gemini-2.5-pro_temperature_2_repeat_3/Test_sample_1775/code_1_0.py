import numpy as np

# --- User-definable parameters ---

# Top gate capacitance per unit area (in Farads per square meter, F/m^2)
Ctg = 1.5e-2 # e.g., 1.5 uF/cm^2 = 1.5e-2 F/m^2

# Back gate capacitance per unit area (in Farads per square meter, F/m^2)
Cbg = 0.5e-2 # e.g., 0.5 uF/cm^2 = 0.5e-2 F/m^2

# Top gate voltage (in Volts, V)
Vtg = 1.5

# Back gate voltage (in Volts, V)
Vbg = 2.0

# --- Calculation ---

# The displacement field (D) is the sum of the contributions from the top and back gates.
# D = C_tg * V_tg + C_bg * V_bg
# The unit of the displacement field is Coulombs per square meter (C/m^2).

D_total = Ctg * Vtg + Cbg * Vbg

# --- Output ---
# The code prints the full equation with the values used, showing the step-by-step calculation.

print("Calculation of the total displacement field (D):")
print("D = C_tg * V_tg + C_bg * V_bg")
print(f"D = ({Ctg:.2e} F/m^2) * ({Vtg:.2f} V) + ({Cbg:.2e} F/m^2) * ({Vbg:.2f} V)")
print(f"D = {Ctg * Vtg:.4f} C/m^2 + {Cbg * Vbg:.4f} C/m^2")
print(f"Total D = {D_total:.4f} C/m^2")

# The final result in a simplified format for completion
# In this example case, 1.5e-2 * 1.5 + 0.5e-2 * 2.0 = 0.0225 + 0.01 = 0.0325
# So the answer is 0.0325