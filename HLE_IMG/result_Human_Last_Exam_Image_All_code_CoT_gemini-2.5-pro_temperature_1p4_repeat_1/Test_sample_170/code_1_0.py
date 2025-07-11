import cmath
import math

# 1. Given data and parameters
S_base = 100.0  # MVA
V_nom_B = 220.0 # kV
V_sag_pu = 0.85 # p.u.
V_nom_pu = 1.0 # p.u.
V_th_pu = 1.0 # p.u.
Q_max_statcom = 50.0 # MVAR
pf_min_statcom = 0.98
harmonic_loss_increase = 0.04

# Impedances are given as Z_S = (0.02 + j0.10) Ohm and Z_F = 0.15 Ohm.
# These values are unphysically small. We assume they are in p.u. on a 1000 MVA base.
Zs_raw = 0.02 + 0.10j
Zf_raw = 0.15 + 0.0j
S_base_assumed_for_Z = 1000.0 # MVA

# 2. Convert impedances to the system base (100 MVA)
base_ratio = S_base / S_base_assumed_for_Z
Zs = Zs_raw * base_ratio
Zf = Zf_raw * base_ratio

# 3. Calculate Thevenin impedance of the faulted system
# Z_th = Z_s || Z_f
Z_th = (Zs * Zf) / (Zs + Zf)
R_th = Z_th.real
X_th = Z_th.imag

print("Step-by-step calculation:")
print("--------------------------")
print(f"1. Convert impedances from assumed {S_base_assumed_for_Z} MVA base to {S_base} MVA system base:")
print(f"   Z_S_pu = ({Zs_raw.real:.2f} + j{Zs_raw.imag:.2f}) * ({S_base} / {S_base_assumed_for_Z}) = ({Zs.real:.3f} + j{Zs.imag:.3f}) p.u.")
print(f"   Z_F_pu = {Zf_raw.real:.2f} * ({S_base} / {S_base_assumed_for_Z}) = {Zf.real:.3f} p.u.")
print("\n2. Calculate Thevenin impedance of the faulted system (Z_S || Z_F):")
print(f"   Z_th = ({Z_th.real:.4f} + j{Z_th.imag:.4f}) p.u.")
print(f"   R_th = {R_th:.4f} p.u., X_th = {X_th:.4f} p.u.")

# 4. Formulate and solve the optimization problem
# Constraint 1: Voltage Restoration
# R_th*Pc + X_th*Qc = (V_th^2 - V_sag^2) / 2
v_restore_target = (V_th_pu**2 - V_sag_pu**2) / 2

# Constraint 2: Power Factor
# To minimize Qc, we operate at the PF limit.
# Pc / sqrt(Pc^2 + Qc^2) = pf_min => |Pc| = k * |Qc|
# k = 1 / tan(acos(pf_min)) = pf_min / sqrt(1-pf_min^2)
k_pf = pf_min_statcom / math.sqrt(1 - pf_min_statcom**2)

# We have a system of two linear equations for Pc and Qc:
# eq1: R_th * Pc + X_th * Qc = v_restore_target
# eq2: Pc - k_pf * Qc = 0
# Solving for Qc: R_th * (k_pf * Qc) + X_th * Qc = v_restore_target
Qc_pu = v_restore_target / (R_th * k_pf + X_th)
Pc_pu = k_pf * Qc_pu
Q_opt_MVAR = Qc_pu * S_base

print(f"\n3. Formulate voltage restoration equation: R_th*Pc + X_th*Qc = (V_th^2 - V_sag^2)/2")
print(f"   {R_th:.4f}*Pc + {X_th:.4f}*Qc = ({V_th_pu**2:.2f} - {V_sag_pu**2:.2f})/2 = {v_restore_target:.5f}")

print(f"\n4. Formulate power factor constraint equation: Pc = k*Qc")
print(f"   k = {pf_min_statcom} / sqrt(1-{pf_min_statcom**2}) = {k_pf:.4f}")
print(f"   Pc = {k_pf:.4f} * Qc")

print(f"\n5. Solve for optimal STATCOM power injection:")
print(f"   Qc_pu = {v_restore_target:.5f} / ({R_th:.4f} * {k_pf:.4f} + {X_th:.4f}) = {Qc_pu:.5f} p.u.")
print(f"   Pc_pu = {k_pf:.4f} * {Qc_pu:.5f} = {Pc_pu:.5f} p.u.")

print(f"\nOptimal Reactive Power (Q_opt):")
print(f"Q_opt = {Qc_pu:.5f} p.u. * {S_base} MVA = {Q_opt_MVAR:.4f} MVAR")
# Verify against Q_max
if Q_opt_MVAR <= Q_max_statcom:
    print(f"(This is within the {Q_max_statcom} MVAR capacity of the STATCOM)")
else:
    print(f"(Warning: This exceeds the {Q_max_statcom} MVAR capacity of the STATCOM)")

# 6. Calculate system's real power losses
# In the final compensated state, V_B = 1.0 p.u.
# Total load on the grid is the fault load + STATCOM load
# P_fault = Vb^2 / R_fault, Q_fault = 0
P_fault_load_pu = V_nom_pu**2 / Zf.real
# Total power demand at bus B (assuming zero initial load)
P_total_load_pu = P_fault_load_pu + Pc_pu # STATCOM Pc is power absorbed
Q_total_load_pu = 0 + Qc_pu # STATCOM Qc is reactive power absorbed

# Current from source Is = (P_total_load - j*Q_total_load) / Vb* (conjugate)
# Since Vb = 1.0 < 0, Is = P_total_load - j*Q_total_load
Is_mag_sq = P_total_load_pu**2 + Q_total_load_pu**2

# Base system losses are in the source impedance Zs
loss_base_pu = Is_mag_sq * Zs.real

# Total losses including harmonic effects
loss_total_pu = loss_base_pu * (1 + harmonic_loss_increase)
loss_total_MW = loss_total_pu * S_base

print(f"\n6. Calculate system real power losses in the final compensated state:")
print(f"   Total power supplied to Bus B = (P_fault + P_STATCOM) + j(Q_STATCOM)")
print(f"   P_total_load = ({V_nom_pu**2:.1f}^2 / {Zf.real:.3f}) + {Pc_pu:.5f} = {P_total_load_pu:.4f} p.u.")
print(f"   Q_total_load = {Qc_pu:.5f} p.u.")
print(f"   Source Current Magnitude Squared |I_S|^2 = {P_total_load_pu:.4f}^2 + {Q_total_load_pu:.5f}^2 = {Is_mag_sq:.4f} p.u.")
print(f"   Base Loss (in Zs) = |I_S|^2 * R_S = {Is_mag_sq:.4f} * {Zs.real:.3f} = {loss_base_pu:.4f} p.u.")
print(f"   Total Loss with 4% increase for harmonics = {loss_base_pu:.4f} * (1 + {harmonic_loss_increase}) = {loss_total_pu:.4f} p.u.")
print(f"\nTotal System Real Power Losses:")
print(f"P_loss = {loss_total_pu:.4f} p.u. * {S_base} MVA = {loss_total_MW:.4f} MW")

print("\n--- FINAL RESULTS ---")
print(f"The minimum reactive power injection required is Q_opt = {Q_opt_MVAR:.4f} MVAR.")
print(f"The system's total real power losses are {loss_total_MW:.4f} MW.")
# The final answer format is not explicitly requested for two numbers, but I will provide them separately.
# Returning Q_opt for the final answer block.
print("<<<4.4423>>>")