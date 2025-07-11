import math

# Step 1: Define system parameters and assumptions
S_base = 100.0  # MVA
V_B_nom = 220.0  # kV
V_A = 1.0  # p.u. (Voltage at external grid, Bus A)

# Per-unit values based on assumptions
# Assumed corrected transmission line impedance
R = 0.2  # p.u.
X = 1.0  # p.u.
Z_S_mag = math.sqrt(R**2 + X**2)

# Fault impedance
Z_F = 0.15  # p.u.

# STATCOM constraints
Q_max_MVAR = 50.0
Q_max_pu = Q_max_MVAR / S_base  # 0.5 p.u.
PF_min = 0.98
# tan(phi) = Q/P, PF = cos(phi). |P/Q| <= tan(acos(PF_min))
P_Q_ratio_max = math.tan(math.acos(PF_min))

# Harmonic loss factor
k_h = 1.04

# Step 2: Analyze pre-compensation state to find P_net
V_B_fault = 0.85  # p.u.

# Fault power at V_B_fault
P_fault_pre = V_B_fault**2 / Z_F

# Using the voltage drop equation: |V_A|-|V_B| = (P*R + Q*X)/|V_B|
# We assume Q_net = 0 (unity power factor for the net load/generation)
# The term (P_net*R + Q_net*X) can be isolated.
# P = P_net + P_fault_pre, Q = Q_net = 0
# |V_A| - V_B_fault = ((P_net + P_fault_pre)*R + 0*X) / V_B_fault
# P_net*R = (|V_A| - V_B_fault) * V_B_fault - P_fault_pre * R
P_net = ((V_A - V_B_fault) * V_B_fault - P_fault_pre * R) / R
# A negative P_net means net power injection from the wind farm.

# Step 3 & 4: Solve for optimal power injection (Q_opt)
V_B_target = 1.0  # p.u.

# The post-compensation voltage equation is:
# |V_A| - V_B_target = ((P_net + P_fault_post - P_comp)*R + (Q_net - Q_comp)*X) / V_B_target
# P_fault_post = V_B_target^2 / Z_F
P_fault_post = V_B_target**2 / Z_F
# Q_net = 0

# With V_A = V_B_target = 1.0, the left side is 0.
# 0 = (P_net + P_fault_post - P_comp)*R + (0 - Q_comp)*X
# 0 = P_net*R + P_fault_post*R - P_comp*R - Q_comp*X
# P_comp*R + Q_comp*X = P_net*R + P_fault_post*R

# We already found P_net*R.
# P_comp*R + Q_comp*X = ((V_A - V_B_fault) * V_B_fault - P_fault_pre * R) + P_fault_post*R
# This simplifies to: Q_comp = C - (R/X)*P_comp, where C is a constant.
constant_term = (V_A - V_B_fault) * V_B_fault * (1/X) + (P_fault_post - P_fault_pre) * (R/X)

# We have Q_comp = constant_term - (R/X)*P_comp
# To minimize Q_comp, we must maximize P_comp.
# Constraint: P_comp <= P_Q_ratio_max * Q_comp
# P_comp <= P_Q_ratio_max * (constant_term - (R/X)*P_comp)
# P_comp * (1 + P_Q_ratio_max * R/X) <= P_Q_ratio_max * constant_term
P_opt = (P_Q_ratio_max * constant_term) / (1 + P_Q_ratio_max * R / X)

# Now find Q_opt using the relationship
Q_opt = constant_term - (R / X) * P_opt

# Step 5: Calculate system losses
# Total power drawn from the grid at Bus B in the compensated state
P_B_total = P_net + P_fault_post - P_opt
Q_B_total = 0 - Q_opt # Q_net is 0

# Squared magnitude of current in the line
S_B_total_sq = P_B_total**2 + Q_B_total**2
I_sq = S_B_total_sq / V_B_target**2

# Line losses without harmonics
P_loss_line = I_sq * R

# Total losses with 4% increase from harmonics
P_loss_total = P_loss_line * k_h

# Convert final values from p.u. to physical units
Q_opt_MVAR = Q_opt * S_base
P_loss_total_MW = P_loss_total * S_base

# --- Output the results ---
print("--- System Parameters and Assumptions ---")
print(f"Assumed Line Impedance Z_S = {R} + j{X} p.u.")
print(f"Fault Impedance Z_F = {Z_F} p.u.")
print(f"Base Power S_base = {S_base} MVA")
print("\n--- Intermediate Calculations ---")
print(f"Net Active Power from Wind Farm (P_net) = {-P_net:.3f} p.u. = {-P_net*S_base:.1f} MW")
print(f"STATCOM Active Power Injection (P_opt) = {P_opt:.3f} p.u. = {P_opt*S_base:.1f} MW")
print("\n--- Final Results ---")
print("Optimization Objective: Minimize reactive power injection from STATCOM.")
print("Solution:")
print(f"The minimum reactive power required from the MMCC STATCOM is:")
print(f"Q_opt = {Q_opt:.4f} p.u. = {Q_opt_MVAR:.2f} MVAR")
print("\nThe system's real power losses under these conditions are:")
print(f"P_loss = {P_loss_total:.4f} p.u. = {P_loss_total_MW:.2f} MW")

# Print the final equations with numbers for clarity
print("\n--- Calculation Summary ---")
print(f"1. Required STATCOM power relation: Q_comp = {constant_term:.4f} - ({R/X:.2f}) * P_comp")
print(f"2. Optimal reactive power: Q_opt = {Q_opt:.4f} p.u.")
print(f"3. Total system loss: P_loss = |I|^2 * R * k_h = {I_sq:.4f} * {R:.2f} * {k_h:.2f} = {P_loss_total:.4f} p.u.")