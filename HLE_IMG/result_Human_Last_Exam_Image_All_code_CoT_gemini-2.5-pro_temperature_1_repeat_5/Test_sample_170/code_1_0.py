import math

# Step 1: Define system parameters in per-unit (pu)
S_base_MVA = 100.0  # MVA
V_base_kV = 220.0   # kV
Q_max_MVAR = 50.0   # MVAR
V_fault_pu = 0.85   # Per-unit voltage during fault
V_nominal_pu = 1.0  # Per-unit nominal voltage

# System impedance (assuming given values are in pu)
Z_S_pu = 0.02 + 0.10j
R_S_pu = Z_S_pu.real
X_S_pu = Z_S_pu.imag

# STATCOM constraints
PF_statcom = 0.98
harmonic_loss_increase = 0.04

# Step 2: Model the voltage drop
# The voltage drop from nominal (1.0 pu) to the fault level (0.85 pu)
delta_V_pu = V_nominal_pu - V_fault_pu

# The voltage drop is caused by the effective load (P_L_eff + jQ_L_eff)
# delta_V = R_S * P_L_eff + X_S * Q_L_eff
# 0.15 = 0.02 * P_L_eff + 0.10 * Q_L_eff

# Step 3: Make an assumption to find the effective load
# Assume the effective load power factor angle is the same as the impedance angle
impedance_angle = math.atan(X_S_pu / R_S_pu)
# So, Q_L_eff / P_L_eff = tan(impedance_angle) = X_S / R_S
# Q_L_eff = P_L_eff * (X_S_pu / R_S_pu)
# Substituting this into the voltage drop equation:
# delta_V = R_S * P_L_eff + X_S * (P_L_eff * X_S / R_S) = P_L_eff * (R_S + X_S^2 / R_S)
# P_L_eff = delta_V / (R_S + X_S**2 / R_S)
P_L_eff_pu = delta_V_pu / (R_S_pu + X_S_pu**2 / R_S_pu)
Q_L_eff_pu = P_L_eff_pu * (X_S_pu / R_S_pu)

# Step 4: Formulate and solve the optimization problem for the STATCOM
# Objective: Minimize Q_comp
# Constraint (Voltage Restoration): R_S*P_comp + X_S*Q_comp = delta_V
# Constraint (Power Factor): We interpret PF > 0.98 for a STATCOM as |Q_comp|/|S_comp| > 0.98,
# which means the device is almost purely reactive.
# |P_comp| / |Q_comp| < tan(acos(|Q|/|S|)) = tan(acos(0.98))
max_P_over_Q_ratio = math.tan(math.acos(PF_statcom))

# To minimize Q_comp, we must maximize P_comp, so we use the limit:
# P_comp = max_P_over_Q_ratio * Q_comp
# Substitute into the voltage restoration constraint:
# R_S * (max_P_over_Q_ratio * Q_comp) + X_S * Q_comp = delta_V
# Q_comp * (R_S * max_P_over_Q_ratio + X_S) = delta_V
Q_opt_pu = delta_V_pu / (R_S_pu * max_P_over_Q_ratio + X_S_pu)
P_comp_pu = max_P_over_Q_ratio * Q_opt_pu

# Convert optimal values from pu to MW and MVAR
Q_opt_MVAR = Q_opt_pu * S_base_MVA
P_comp_MW = P_comp_pu * S_base_MVA

# Step 5: Calculate system losses in the final compensated state
# Final net power drawn at Bus B from the grid
P_B_final_pu = P_L_eff_pu - P_comp_pu
Q_B_final_pu = Q_L_eff_pu - Q_opt_pu

# Final line current squared (|I|^2 = |S|^2 / |V|^2)
# Since V_B is restored to 1.0 pu, |V|^2 = 1.0
I_sq_pu = P_B_final_pu**2 + Q_B_final_pu**2

# Line losses (I^2 * R)
P_loss_line_pu = I_sq_pu * R_S_pu

# Total system losses before harmonic effects
# These are the line losses plus the STATCOM's own real power consumption
P_loss_system_pu = P_loss_line_pu + P_comp_pu

# Account for the 4% increase due to harmonics
P_loss_total_pu = P_loss_system_pu * (1 + harmonic_loss_increase)

# Convert total losses to MW
P_loss_total_MW = P_loss_total_pu * S_base_MVA

# Step 6: Print the final results
print("--- Optimization Results ---")
print(f"The minimum required reactive power from the STATCOM is:")
print(f"Q_opt = {Q_opt_MVAR:.2f} MVAR")
print("\nNote: This required value exceeds the STATCOM's maximum capacity of 50 MVAR.")
print("Therefore, restoring the voltage to the nominal value is not feasible with the given equipment.\n")
print("--- System Power Loss Calculation ---")
print("Assuming the voltage is restored as per the optimization goal, the system losses are:")
print(f"STATCOM Real Power Consumption: P_comp = {P_comp_MW:.2f} MW")
print(f"Total System Real Power Losses (including harmonic effects): P_loss = {P_loss_total_MW:.2f} MW")

print("\n--- Final Equation for Losses ---")
print("P_loss_total = ( |(P_L_eff - P_comp) + j(Q_L_eff - Q_opt)|^2 * R_s + P_comp ) * (1 + Harmonic_Factor)")
print(f"P_loss_total = ( |({P_L_eff_pu:.4f} - {P_comp_pu:.4f}) + j({Q_L_eff_pu:.4f} - {Q_opt_pu:.4f})|^2 * {R_S_pu:.2f} + {P_comp_pu:.4f} ) * (1 + {harmonic_loss_increase:.2f})")
print(f"P_loss_total = ( |({P_B_final_pu:.4f}) + j({Q_B_final_pu:.4f})|^2 * {R_S_pu:.2f} + {P_comp_pu:.4f} ) * {1+harmonic_loss_increase:.2f}")
print(f"P_loss_total = ( {I_sq_pu:.6f} * {R_S_pu:.2f} + {P_comp_pu:.4f} ) * {1+harmonic_loss_increase:.2f}")
print(f"P_loss_total = ( {P_loss_line_pu:.6f} + {P_comp_pu:.4f} ) * {1+harmonic_loss_increase:.2f}")
print(f"P_loss_total = {P_loss_system_pu:.4f} pu * {1+harmonic_loss_increase:.2f} = {P_loss_total_pu:.4f} pu")
print(f"P_loss_total_MW = {P_loss_total_pu:.4f} * {S_base_MVA:.0f} MVA = {P_loss_total_MW:.2f} MW")