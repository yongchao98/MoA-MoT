import math

# Step 1: Define the given parameters from the problem description.
# The voltage at Bus 4 is to be maintained at 400 kV.
V_bus4_kV = 400.0  # Voltage in kilovolts (kV)

# The reactance of the shunt capacitor bank is 5 Ohms.
Xc_ohms = 5.0  # Reactance in Ohms (Ω)

# Step 2: Calculate the reactive power supplied by the capacitor bank.
# The formula for reactive power (Qc) from a shunt capacitor is Qc = V^2 / Xc.
# When V is in kV and Xc is in Ohms, the resulting Qc is in MVAR (Mega-VARs).
# This is because (kV^2 / Ω) = ((10^3 V)^2 / Ω) = (10^6 V^2 / Ω) = 10^6 VAR = 1 MVAR.
Qc_MVAR = (V_bus4_kV ** 2) / Xc_ohms

# Step 3: Print the analysis and the result.
# The reactive power needed to maintain the voltage is equal to the
# reactive power supplied by the specified capacitor bank at that voltage.
print("Reactive Power Compensation Analysis at Bus 4")
print("---------------------------------------------")
print("The goal is to find the reactive power (Qc) needed to maintain the voltage at Bus 4 at 400 kV.")
print("This compensation is provided by a shunt capacitor with a given reactance.")
print("\nFormula for reactive power supplied by a capacitor:")
print("Qc (MVAR) = (Voltage (kV))^2 / Reactance (Ω)")
print("\nGiven values:")
print(f"Voltage at Bus 4 (V): {V_bus4_kV} kV")
print(f"Capacitor Reactance (Xc): {Xc_ohms} Ω")
print("\nCalculation:")
print(f"Qc = ({V_bus4_kV})^2 / {Xc_ohms}")
print(f"Qc = {V_bus4_kV**2} / {Xc_ohms}")
print(f"Qc = {Qc_MVAR:.0f} MVAR")

print("\nTherefore, the reactive power needed to maintain the bus voltage at 400 kV is 32000 MVAR.")

# Final answer block
final_answer = int(Qc_MVAR)
print(f"\n<<<{final_answer}>>>")