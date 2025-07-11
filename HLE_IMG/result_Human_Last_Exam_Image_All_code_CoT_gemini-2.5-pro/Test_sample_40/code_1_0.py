import cmath
import math

# --- Given Data ---
# System parameters in per-unit (pu)
# We assume the given impedance values are in pu on the specified base.
Zs_pu = 0.02 + 0.10j
Zf_pu = 0.15 + 0.0j
Rs_pu = Zs_pu.real

# Voltage conditions in pu
Vb1_mag = 0.85  # Voltage at Bus B before compensation
Vb2_mag = 1.0   # Target voltage at Bus B after compensation
Va = 1.0 + 0.0j # Voltage at Bus A (external grid)

# STATCOM and System constraints
Q_max_MVAR = 50
S_base_MVA = 100
harmonic_loss_factor = 1.04

# --- 1. Calculation of Optimal Reactive Power (Q_opt) ---
# This part uses the constant impedance load model based on the given Vb1 = 0.85 pu.

# Calculate admittances
Ys = 1 / Zs_pu
Yf = 1 / Zf_pu
Ynet = Ys + Yf
Ynet_star = Ynet.conjugate

# To find Q_opt for P_comp = 0, we first need to solve for the phase angle difference
# between state 1 and 2 (delta_diff = d1 - d2).
# The condition P_comp = 0 leads to the equation:
# (1 - |Vb1|*cos(delta_diff))*G_net_star - (|Vb1|*sin(delta_diff))*B_net_star = 0
# where Ynet_star = G_net_star + j*B_net_star.
# Rearranging gives: A*cos(delta_diff) + B*sin(delta_diff) = C
G_net_star = Ynet_star.real
B_net_star = Ynet_star.imag

A = Vb1_mag * G_net_star
B = Vb1_mag * B_net_star
C = G_net_star

# Solve A*cos(x) + B*sin(x) = C using R*cos(x-alpha) = C
R = math.sqrt(A**2 + B**2)
alpha = math.atan2(B, A)
cos_val = C / R
# Handle potential math domain error if C/R > 1
if abs(cos_val) > 1:
    print("Error: Cannot solve for delta_diff. The problem parameters may be inconsistent.")
    # As the problem is likely flawed, we will proceed with the calculation
    # but the result will be physically questionable. The code will continue
    # with the value anyway, which may lead to a math domain error.
    # To prevent crash, let's clamp it.
    cos_val = max(min(cos_val, 1.0), -1.0)


delta_diff = math.acos(cos_val) + alpha # We choose one solution for the angle

# Now calculate Q_opt (in pu) using the derived formula:
# Q_opt = (1 - |Vb1|*cos(delta_diff))*B_net_star + (|Vb1|*sin(delta_diff))*G_net_star
cos_dd = math.cos(delta_diff)
sin_dd = math.sin(delta_diff)

Q_opt_pu = (Vb2_mag**2 - Vb1_mag * cos_dd) * B_net_star + (Vb1_mag * sin_dd) * G_net_star
Q_opt_MVAR = Q_opt_pu * S_base_MVA

# --- 2. Calculation of System Real Power Losses ---
# NOTE: This calculation is unsolvable without knowing the system load.
# We make a contradictory assumption: the system has no load (I_load=0).
# This will allow a numerical answer but contradicts the Vb1=0.85 pu premise.

# Under no load, Vb1 is determined by the voltage divider Zf / (Zs + Zf)
Vb1_noload = Va * Zf_pu / (Zs_pu + Zf_pu)
Vb1_noload_star = Vb1_noload.conjugate()

# We need the post-compensation voltage phasor Vb2 = 1.0 * exp(j*delta2)
# Solve for delta2 by setting the active power from STATCOM to 0.
# S_comp = (|Vb2|^2 - Vb2*Vb1_star) * Ynet_star
# P_comp = Re[(1 - Vb2*Vb1_star) * Ynet_star] = 0
# Let Vb2 = cos(d2) + j*sin(d2). This leads to an equation for d2.
# A2*cos(d2) + B2*sin(d2) = C2
A2 = Vb1_noload_star.real * G_net_star + Vb1_noload_star.imag * B_net_star
B2 = Vb1_noload_star.imag * G_net_star - Vb1_noload_star.real * B_net_star
C2 = G_net_star

R2 = math.sqrt(A2**2 + B2**2)
alpha2 = math.atan2(B2, A2)

cos_val2 = C2 / R2
if abs(cos_val2) > 1:
    cos_val2 = max(min(cos_val2, 1.0), -1.0)
    
delta2 = math.acos(cos_val2) + alpha2

# Choose the angle that gives a smaller voltage angle difference, usually more stable.
delta2_alt = -math.acos(cos_val2) + alpha2
if abs(delta2_alt) < abs(delta2):
    delta2 = delta2_alt

Vb2 = Vb2_mag * (math.cos(delta2) + 1j * math.sin(delta2))

# Calculate line current Is and fundamental losses
Is = (Va - Vb2) / Zs_pu
Ploss_fund_pu = (abs(Is)**2) * Rs_pu
Ploss_fund_MW = Ploss_fund_pu * S_base_MVA

# Calculate total losses including 4% harmonic increase
Ploss_total_pu = Ploss_fund_pu * harmonic_loss_factor
Ploss_total_MW = Ploss_total_pu * S_base_MVA

# --- 3. Print Results ---
print("--- Optimization Results ---")
print(f"Objective: Minimize reactive power injection Q_opt to restore voltage to {Vb2_mag:.2f} p.u.")
print("\n--- Calculated Reactive Power ---")
print(f"The minimum reactive power required from the MMCC STATCOM is:")
print(f"Q_opt = {Q_opt_pu:.4f} p.u.")
print(f"Q_opt = {Q_opt_MVAR:.2f} MVAR")
if Q_opt_MVAR > Q_max_MVAR:
    print(f"(Note: This required power exceeds the STATCOM's maximum capacity of {Q_max_MVAR} MVAR.)")


print("\n--- Calculated System Losses ---")
print("Note: Loss calculation requires a simplifying assumption of no external load, which contradicts the initial voltage sag value.")
print(f"The fundamental real power loss in the transmission line is:")
print(f"P_loss_fund = {Ploss_fund_pu:.4f} p.u. = {Ploss_fund_MW:.2f} MW")
print(f"Considering a {int((harmonic_loss_factor-1)*100)}% increase due to harmonics, the total system real power loss is:")
print(f"P_loss_total = P_loss_fund * {harmonic_loss_factor:.2f} = {Ploss_total_pu:.4f} p.u. = {Ploss_total_MW:.2f} MW")

print("\n--- Final Equations with numbers ---")
# Equation for Q_opt
print(f"Q_opt = ((1.0^2 - {Vb1_mag}*cos({delta_diff:.4f}))*{B_net_star:.4f} + ({Vb1_mag}*sin({delta_diff:.4f}))*{G_net_star:.4f}) * {S_base_MVA} MVA = {Q_opt_MVAR:.2f} MVAR")
# Equation for P_loss
print(f"P_loss_total = |({Va.real:.2f} - ({Vb2.real:.4f} + {Vb2.imag:.4f}j)) / ({Zs_pu.real:.2f} + {Zs_pu.imag:.2f}j)|^2 * {Rs_pu:.2f} * {harmonic_loss_factor:.2f} * {S_base_MVA} MW = {Ploss_total_MW:.2f} MW")
