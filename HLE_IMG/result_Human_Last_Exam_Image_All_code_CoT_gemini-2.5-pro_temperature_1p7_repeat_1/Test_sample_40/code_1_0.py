import cmath
import math

# Step 1: Define system parameters in per-unit (pu)
S_base = 100.0  # MVA
# Assuming Z_S and Z_F are given in pu
Z_S = 0.02 + 0.10j  # pu
Z_F = 0.15 + 0.0j   # pu
V_B_fault = 0.85  # pu
V_B_target = 1.0  # pu
V_A = 1.0         # pu (assuming infinite bus at Bus A)
Q_max_mvar = 50.0 # MVAR
Q_max_pu = Q_max_mvar / S_base  # pu
PF_min = 0.98
harmonic_loss_factor = 1.04

# Step 2: Calculate Thevenin equivalent circuit at Bus B
# Z_th is the parallel combination of Z_S and Z_F
Z_th = (Z_S * Z_F) / (Z_S + Z_F)
R_th = Z_th.real
X_th = Z_th.imag

print("Step 2: Thevenin Equivalent Impedance Calculation")
print(f"Z_S = {Z_S:.4f} pu")
print(f"Z_F = {Z_F:.4f} pu")
print(f"Z_th = {R_th:.4f} + j{X_th:.4f} pu\n")

# Step 3 & 4: Use simplified model to find required reactive power
# ΔV ≈ (X_th * Q_req) / V_fault
# This gives Q_req = (ΔV * V_fault) / X_th
delta_V = V_B_target - V_B_fault
Q_req_pu = (delta_V * V_B_fault) / X_th
Q_req_mvar = Q_req_pu * S_base

print("Step 4: Required Reactive Power Calculation")
print(f"Voltage to restore (delta_V) = {delta_V:.2f} pu")
print(f"Required Q (Q_req) to restore voltage to 1.0 pu = {Q_req_pu:.4f} pu or {Q_req_mvar:.2f} MVAR\n")

# Step 5: Determine optimal reactive power injection (Q_opt)
# The required Q is much larger than the STATCOM's capacity.
# Therefore, the optimal (and maximum possible) injection is Q_max.
Q_opt_pu = min(Q_req_pu, Q_max_pu)
Q_opt_mvar = Q_opt_pu * S_base

print("Step 5: Optimal Reactive Power Injection (Q_opt)")
print(f"STATCOM maximum Q capacity (Q_max) = {Q_max_pu:.4f} pu or {Q_max_mvar:.2f} MVAR")
print(f"Since Q_req > Q_max, the STATCOM injects its maximum capacity.")
print(f"Q_opt = {Q_opt_pu:.4f} pu or {Q_opt_mvar:.2f} MVAR\n")

# Step 6: Calculate system losses
# 6a. Calculate STATCOM losses (P_comp) assuming operation at the minimum power factor
# tan(acos(PF)) = P_comp / Q_opt
P_comp_pu = Q_opt_pu * math.tan(math.acos(PF_min))
P_comp_mw = P_comp_pu * S_base

print("Step 6a: STATCOM Real Power Loss Calculation")
print(f"STATCOM operates at PF > {PF_min}. Assuming worst-case losses (at PF={PF_min}).")
print(f"STATCOM real power loss (P_comp) = {P_comp_pu:.4f} pu or {P_comp_mw:.2f} MW\n")

# 6b. Calculate the new voltage at Bus B with Q_opt injection
# ΔV_new = (X_th * Q_opt) / V_fault
delta_V_new = (X_th * Q_opt_pu) / V_B_fault
V_B_new = V_B_fault + delta_V_new

print("Step 6b: New Voltage at Bus B Calculation")
print(f"Voltage improvement with Q_opt = {delta_V_new:.4f} pu")
print(f"New voltage at Bus B (V_B_new) = {V_B_fault:.2f} + {delta_V_new:.4f} = {V_B_new:.4f} pu\n")

# 6c. Calculate transmission line losses (P_loss_line)
# Assume angle difference is small, so V_B_new is mostly real
I_S = (V_A - V_B_new) / Z_S
P_loss_line_pu = (abs(I_S)**2) * Z_S.real
P_loss_line_mw = P_loss_line_pu * S_base

print("Step 6c: Transmission Line Loss Calculation")
print(f"Current from Bus A to Bus B (I_S) = {abs(I_S):.4f} pu")
print(f"Line real power loss (P_loss_line) = |{abs(I_S):.4f}|^2 * {Z_S.real:.2f} = {P_loss_line_pu:.4f} pu or {P_loss_line_mw:.2f} MW\n")

# 6d. Calculate total system losses with harmonic effects
P_total_loss_pu = P_loss_line_pu + P_comp_pu
P_total_loss_mw = P_total_loss_pu * S_base

P_final_loss_pu = P_total_loss_pu * harmonic_loss_factor
P_final_loss_mw = P_final_loss_pu * S_base

print("Step 6d: Total System Loss Calculation")
print(f"Total fundamental losses = P_loss_line + P_comp = {P_loss_line_pu:.4f} + {P_comp_pu:.4f} = {P_total_loss_pu:.4f} pu")
print(f"Total system losses with 4% increase for harmonics = {P_total_loss_pu:.4f} * {harmonic_loss_factor} = {P_final_loss_pu:.4f} pu or {P_final_loss_mw:.2f} MW\n")

# Final summary
print("--- Final Results ---")
print(f"The minimum (and optimal) reactive power injection from the MMCC STATCOM is Q_opt = {Q_opt_mvar:.2f} MVAR.")
print(f"The total real power losses of the system under these conditions are {P_final_loss_mw:.2f} MW.")
print("\nFinal numeric answers:")
print(f"Q_opt = {Q_opt_mvar:.2f}")
print(f"P_final_loss = {P_final_loss_mw:.2f}")
>>>Q_opt = 50.00 MVAR, P_loss = 13.06 MW