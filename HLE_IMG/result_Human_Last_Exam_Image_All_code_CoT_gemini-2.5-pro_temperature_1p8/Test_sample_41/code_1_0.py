import math

# --- Given Parameters ---
V_ac_nominal = 500e3  # V, Assumed AC bus nominal voltage (line-to-line)
L = 0.1  # H, Assumed AC system Thevenin inductance
voltage_drop_pct = 2.5  # %
P_inverter = 254.97e6  # W, Active power from the inverter
HD_3_pct = 5.0  # %, 3rd harmonic current distortion
HD_5_pct = 10.0 # %, 5th harmonic current distortion

# --- Assumptions ---
f = 60  # Hz, Assumed system frequency

# --- Step 1: Calculate Reactive Power Compensation for Voltage Drop (Q_fault) ---

# The reactive power change (delta_Q) required to produce a small voltage
# change (delta_V) at a bus is approximately delta_Q = (delta_V / V) * (V^2 / X),
# where X is the system reactance and V^2/X is the short-circuit MVA.
# delta_V/V is the per-unit voltage change.

# Calculate Thevenin reactance of the AC system
X_th = 2 * math.pi * f * L

# Calculate the reactive power compensation needed to restore the voltage drop
delta_V_pu = voltage_drop_pct / 100.0
# The formula simplifies to delta_Q = delta_V_pu * (V_ac_nominal^2 / X_th)
q_fault = delta_V_pu * (V_ac_nominal**2) / X_th

# Convert to MVar for readability
q_fault_mvar = q_fault / 1e6

# --- Step 2: Calculate Reactive Power Compensation for Harmonics (Q_harmonics) ---

# This is interpreted as the distortion power (D), which is calculated as D = S1 * THDi
# where S1 is the fundamental apparent power and THDi is the total harmonic current distortion.

# First, estimate the fundamental reactive power (Q1) consumed by the inverter.
# A common rule of thumb for 12-pulse LCC converters is Q1 is 50% of P1.
Q1_inverter = 0.5 * P_inverter

# Calculate the fundamental apparent power (S1) of the inverter
S1_inverter = math.sqrt(P_inverter**2 + Q1_inverter**2)

# Calculate the Total Harmonic Distortion (THD) for current (i)
hd_3 = HD_3_pct / 100.0
hd_5 = HD_5_pct / 100.0
thd_i = math.sqrt(hd_3**2 + hd_5**2)

# Calculate the distortion power (D)
q_harmonics = S1_inverter * thd_i

# Convert to MVar for readability
q_harmonics_mvar = q_harmonics / 1e6

# --- Step 3: Calculate Total Reactive Power Compensation ---

q_total_mvar = q_fault_mvar + q_harmonics_mvar

# --- Output the results ---
print("Calculation Steps and Final Answer:")
print("-----------------------------------")
print(f"1. Compensation for Voltage Drop (Q_fault):")
print(f"   - Assumed AC Voltage (V): {V_ac_nominal/1e3:.0f} kV")
print(f"   - System Reactance (X_th = 2*pi*f*L): {X_th:.2f} Ohms")
print(f"   - Voltage Drop (delta_V): {voltage_drop_pct:.1f}%")
print(f"   - Q_fault = ({voltage_drop_pct/100:.3f}) * ({V_ac_nominal/1e3:.0f} kV)^2 / {X_th:.2f} Ohms = {q_fault_mvar:.2f} MVar")
print("\n2. Compensation for Harmonics (Q_harmonics):")
print(f"   - Inverter Active Power (P1): {P_inverter/1e6:.2f} MW")
print(f"   - Estimated Inverter Reactive Power (Q1): {Q1_inverter/1e6:.2f} MVar")
print(f"   - Inverter Apparent Power (S1): {S1_inverter/1e6:.2f} MVA")
print(f"   - Total Harmonic Distortion (THDi): {thd_i*100:.2f}%")
print(f"   - Q_harmonics = S1 * THDi = {S1_inverter/1e6:.2f} MVA * {thd_i:.4f} = {q_harmonics_mvar:.2f} MVar")
print("\n3. Total Compensation (Q_total):")
print(f"   Q_total = Q_fault + Q_harmonics")
# The final result line as requested
print(f"   Total Reactive Power Compensation = {q_fault_mvar:.2f} MVar + {q_harmonics_mvar:.2f} MVar = {q_total_mvar:.2f} MVar")
