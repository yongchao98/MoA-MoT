import math

# --- Problem Parameters ---
# The target line-to-line voltage to be maintained at Bus 4.
V_L_kV = 400.0  # in kV
# The reactance of the shunt capacitor bank providing the compensation.
Xc = 5.0  # in Ohms

# --- Analysis ---
# The problem asks for the reactive power needed to maintain the voltage at 400 kV,
# and specifies that this compensation is provided by a shunt capacitor of 5 Ohms.
# The reactive power (Q) supplied by a shunt capacitor is given by the formula Q = V^2 / Xc.
# Here, V is the line-to-line voltage.

# Convert the voltage from kilovolts (kV) to volts (V) for SI unit consistency.
V_L_V = V_L_kV * 1000

# Calculate the reactive power in VAR (Volt-Ampere Reactive).
Q_VAR = V_L_V**2 / Xc

# Convert the result to MVAR (MegaVAR) for convenience, as is standard in power systems.
Q_MVAR = Q_VAR / 1_000_000

# --- Output Results ---
print("--- Reactive Power Compensation Analysis ---")
print("The analysis determines the reactive power needed to maintain the voltage at Bus 4 at 400 kV using a specified shunt capacitor.")
print("\nThe reactive power (Q) supplied by the shunt capacitor is calculated using the formula: Q = V^2 / Xc\n")

print(f"Given values:")
print(f"Target Bus Voltage (V) = {V_L_kV} kV = {int(V_L_V)} V")
print(f"Shunt Capacitor Reactance (Xc) = {Xc} Ω\n")

print("Calculation Steps:")
print(f"Q = ({int(V_L_V)})^2 / {int(Xc)}")
print(f"Q = {int(V_L_V**2)} / {int(Xc)}")
print(f"Q = {int(Q_VAR)} VAR")
print(f"Q = {int(Q_MVAR)} MVAR\n")

print(f"Therefore, the reactive power needed to be supplied by the capacitor bank to maintain the voltage at 400 kV is {int(Q_MVAR)} MVAR.")
