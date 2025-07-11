import math

# --- Given Parameters ---
V1 = 0.78  # Start voltage in Volts
V2 = 0.98  # End voltage in Volts
I2 = 0.445 # Change in current in Amperes
R_load = 50.0  # Load resistance in Ohms
margin = 0.20  # 20% startup margin

# --- Step 1 & 2: Calculate Diode's Dynamic Resistance (rd) ---
# The change in voltage (delta_V)
delta_V = V2 - V1
# The change in current (delta_I), assuming the change is from 0 to I2.
delta_I = I2
# The dynamic resistance (rd) is delta_V / delta_I
rd = delta_V / delta_I

# --- Step 3: Apply the Startup Margin ---
# The target design impedance includes the 20% margin
rd_design = rd * (1 + margin)

# --- Step 4: Calculate the Impedance Transformation Ratio ---
# For optimal power transfer, rd_design = Z_ratio * R_load
# So, Z_ratio = rd_design / R_load
impedance_ratio = rd_design / R_load

# --- Final Output ---
print("Calculation Steps:")
print(f"1. The change in voltage (ΔV) = {V2} V - {V1} V = {delta_V:.3f} V")
print(f"2. The change in current (ΔI) is assumed to be {delta_I} A")
print(f"3. The diode's dynamic resistance (rd) = {delta_V:.3f} V / {delta_I} A = {rd:.3f} Ohms")
print(f"4. Applying the {margin*100}% margin: rd_design = {rd:.3f} * (1 + {margin}) = {rd_design:.3f} Ohms")
print(f"5. The load resistance (R_load) is {R_load} Ohms")

print("\nThe final equation for the impedance transformation ratio (Z_ratio) is:")
print(f"Z_ratio = ((V2 - V1) / I2) * (1 + margin) / R_load")
print(f"Z_ratio = (({V2} - {V1}) / {I2}) * (1 + {margin}) / {R_load}")
print(f"\nResult:")
print(f"The required impedance transformation ratio is: {impedance_ratio}")

# The final answer in the required format
final_answer = impedance_ratio