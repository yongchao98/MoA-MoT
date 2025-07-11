import math

# Step 1: Define initial and target voltages
V_nominal_kV = 400.0  # Nominal line-to-line voltage in kV
voltage_drop_percent = 5.0 # Voltage drop at Bus 4

# Calculate the uncompensated voltage at Bus 4
V_uncomp_kV = V_nominal_kV * (1 - voltage_drop_percent / 100.0)

# The target compensated voltage is the nominal voltage
V_comp_kV = V_nominal_kV

# Convert voltages to Volts for calculation
V_comp_V = V_comp_kV * 1000.0
V_uncomp_V = V_uncomp_kV * 1000.0

# Step 2: Calculate the required voltage correction
delta_V_V = V_comp_V - V_uncomp_V

# Step 3: Define the system Thevenin reactance
# As explained, we interpret the given Xc = 5 Ohms as the system's Thevenin reactance at Bus 4.
X_th_ohm = 5.0

# Step 4: Calculate the required reactive power using the approximation formula
# Q_needed (VAR) = (V * delta_V) / X_th
# We use the target compensated voltage in the formula.
Q_needed_VAR = (V_comp_V * delta_V_V) / X_th_ohm

# Convert the result to Megavars (MVAR)
Q_needed_MVAR = Q_needed_VAR / 1e6

# Step 5: Print the results and the calculation
print("Reactive Power Compensation Analysis for Bus 4\n")
print(f"Nominal Voltage (V_nominal): {V_nominal_kV} kV")
print(f"Voltage Drop: {voltage_drop_percent}%")
print(f"Uncompensated Voltage at Bus 4 (V_uncomp): {V_uncomp_kV} kV")
print(f"Target Compensated Voltage at Bus 4 (V_comp): {V_comp_kV} kV\n")
print(f"Required Voltage Increase (delta_V): {delta_V_V / 1000.0} kV")
print(f"System Thevenin Reactance at Bus 4 (X_th): {X_th_ohm} Ohms\n")

print("The required reactive power (Q_needed) is calculated as:")
print("Q_needed = (V_comp * delta_V) / X_th")
print(f"Q_needed = ({int(V_comp_V)} V * {int(delta_V_V)} V) / {X_th_ohm} Ω")
print(f"Q_needed = {int(Q_needed_VAR)} VAR")
print(f"Q_needed = {Q_needed_MVAR} MVAR\n")
print("Therefore, the reactive power needed to maintain the bus voltage at 400 kV is 1600 MVAR.")

<<<1600.0>>>