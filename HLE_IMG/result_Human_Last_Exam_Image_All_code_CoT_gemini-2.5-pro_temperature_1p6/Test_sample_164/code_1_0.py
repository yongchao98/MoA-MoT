import math

# Step 1: Define the given parameters from the problem description.
V_bus4_target_kV = 400  # Target line-to-line voltage at Bus 4 in kV
Xc_ohms = 5             # Reactance of the shunt capacitor bank in Ohms

# Step 2: Convert the voltage from kV to V for calculation.
V_bus4_target_V = V_bus4_target_kV * 1000

# Step 3: Calculate the reactive power supplied by the capacitor.
# The formula for total three-phase reactive power (Qc) is V_line^2 / Xc.
# The result will be in VAR (Volt-Ampere Reactive).
Qc_var = (V_bus4_target_V)**2 / Xc_ohms

# Step 4: Convert the result from VAR to MVAR (Mega-VAR).
Qc_mvar = Qc_var / 1_000_000

# Step 5: Display the detailed analysis and the final result.
print("Reactive Power Compensation Analysis")
print("="*40)
print(f"The goal is to find the reactive power needed to maintain the voltage at Bus 4 at {V_bus4_target_kV} kV.")
print(f"This is provided by a shunt capacitor with a reactance of {Xc_ohms} Ω.")
print("\nThe reactive power (Qc) is calculated with the formula: Qc = V^2 / Xc")
print("\nFinal Equation with values:")
# Using format to make numbers more readable.
print(f"Qc (MVAR) = ({int(V_bus4_target_V):,})^2 V^2 / {Xc_ohms} Ω / 1,000,000")
print(f"Qc (MVAR) = {int(Qc_mvar):,}")

print("\nTherefore, the reactive power supplied by the capacitor bank is {:,} MVAR.".format(int(Qc_mvar)))
<<<32000>>>