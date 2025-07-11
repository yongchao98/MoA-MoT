import math

# --- Given Parameters ---
P_inv_MW = 254.97  # Active power at the inverter in MW
V_dc_kV = 500      # HVDC voltage in kV
L_H = 0.1          # Line inductance in Henry
volt_drop_pct = 2.5 # Voltage drop at Bus 6 in percent
h3_dist_pct = 5     # 3rd harmonic distortion in percent
h5_dist_pct = 10    # 5th harmonic distortion in percent

# --- Assumptions ---
f_Hz = 60.0 # Assumed system frequency in Hz

# --- Calculations ---

# 1. Calculate Reactive Power for Harmonics (Q_harmonics)
# Convert distortion percentages to per-unit values
h3_dist_pu = h3_dist_pct / 100.0
h5_dist_pu = h5_dist_pct / 100.0

# Calculate current Total Harmonic Distortion (I_THD)
I_THD = math.sqrt(h3_dist_pu**2 + h5_dist_pu**2)

# Calculate reactive power due to harmonics (Distortion Power D) in MVar
# D is approximated as P * I_THD
Q_harmonics_MVar = P_inv_MW * I_THD

# 2. Calculate Reactive Power for Fault (Q_fault)
# Convert voltage drop percentage to a per-unit value
volt_drop_pu = volt_drop_pct / 100.0

# Convert V_dc from kV to V for calculation
V_V = V_dc_kV * 1000

# Calculate angular frequency (omega)
omega_rad_s = 2 * math.pi * f_Hz

# Calculate the system reactance X_sys from L
X_sys_ohm = omega_rad_s * L_H

# Calculate the reactive power compensation needed to restore the voltage
# Formula: Q_fault = (delta_V/V) * V^2 / X_sys
# We use volt_drop_pu for (delta_V/V) and assume V = V_dc
Q_fault_Var = volt_drop_pu * (V_V**2) / X_sys_ohm

# Convert from Var to MVar
Q_fault_MVar = Q_fault_Var / 1e6

# 3. Calculate Total Reactive Power Compensation
Q_total_MVar = Q_fault_MVar + Q_harmonics_MVar

# --- Output the Results ---
print("Calculation of Total Reactive Power Compensation:")
print("-" * 50)
print(f"1. Reactive Power for Harmonics (Q_harmonics):")
print(f"   Based on the 3rd harmonic distortion of {h3_dist_pct}% and 5th harmonic distortion of {h5_dist_pct}%,")
print(f"   the harmonic reactive power required is {Q_harmonics_MVar:.2f} MVar.")
print("\n2. Reactive Power for Fault-induced Voltage Drop (Q_fault):")
print(f"   To correct a {volt_drop_pct}% voltage drop, the reactive power required is {Q_fault_MVar:.2f} MVar.")
print("-" * 50)
print("Final Equation for Total Compensation:")
print(f"Q_total = Q_fault + Q_harmonics")
print(f"Q_total = {Q_fault_MVar:.2f} MVar + {Q_harmonics_MVar:.2f} MVar = {Q_total_MVar:.2f} MVar")

print(f"\nTherefore, the total reactive power compensation required is {Q_total_MVar:.2f} MVar.")
<<<194.30>>>