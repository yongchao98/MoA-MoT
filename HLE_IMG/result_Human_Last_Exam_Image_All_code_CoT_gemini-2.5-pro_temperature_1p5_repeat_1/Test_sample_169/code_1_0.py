import cmath
import math

# Step 1: Define given constants and problem parameters
S_base = 100e6  # Base Power in VA (100 MVA)
V_nom_target = 0.575e3 # Nominal and target voltage in V (0.575 kV)
Z_WF_ohm = 0.01 + 0.05j # Transmission impedance in Ohms
R_F_ohm = 0.1 # Fault resistance in Ohms
Q_max_mvar = 10.0 # Maximum reactive power injection in MVAR
min_pf = 0.95 # Minimum lagging power factor
harmonic_loss_increase = 0.06 # 6% increase in total system losses

# Step 2: Establish the Per-Unit System
V_base = V_nom_target
Z_base = V_base**2 / S_base

# Step 3: Convert all values to Per-Unit (p.u.)
Z_WF_pu = Z_WF_ohm / Z_base
R_F_pu = R_F_ohm / Z_base
V_W_target_pu = V_nom_target / V_base # This will be 1.0 p.u.
Q_max_pu = Q_max_mvar * 1e6 / S_base

# Step 4: Model the system during the fault
# Total impedance from Bus-W to ground through the fault
Z_total_pu = Z_WF_pu + R_F_pu
R_total_pu = Z_total_pu.real
X_total_pu = Z_total_pu.imag

# Step 5: Calculate required power draw to meet the voltage target
# S_W = |V_W|^2 / Z_total*. Here |V_W| = 1.0 p.u.
S_W_pu = (V_W_target_pu**2) / Z_total_pu.conjugate()
P_W_pu = S_W_pu.real
Q_W_pu = S_W_pu.imag

# Step 6: Account for harmonic losses
# The real power generated must cover the fundamental losses (P_W_pu) plus harmonic losses
P_gen_pu = P_W_pu * (1 + harmonic_loss_increase)

# Step 7: Apply constraints to find the optimal Q_comp
# To minimize Q_comp, we must maximize Q_gen from the generator.
# The maximum Q_gen is limited by the power factor constraint.
# PF = cos(phi) >= 0.95 => tan(phi) <= tan(acos(0.95))
# tan(phi) = Q_gen / P_gen
max_tan_phi = math.tan(math.acos(min_pf))
Q_gen_max_pu = P_gen_pu * max_tan_phi

# From power balance at Bus-W: Q_gen + Q_comp = Q_W
# Therefore, Q_comp = Q_W - Q_gen
# To minimize Q_comp, we use the maximum possible Q_gen
Q_opt_pu = Q_W_pu - Q_gen_max_pu

# Convert the final optimal reactive power back to MVAR
Q_opt_MVAR = Q_opt_pu * (S_base / 1e6)

# Step 8: Print the results and the final equation components
print("--- System Parameters and Formulation ---")
print(f"Base Power (S_base): {S_base/1e6} MVA")
print(f"Base Voltage (V_base): {V_base/1e3} kV")
print(f"Base Impedance (Z_base): {Z_base:.6f} Ohms")
print(f"Target Bus-W Voltage: {V_W_target_pu:.3f} p.u.")
print("\n--- Power Calculation to meet Voltage Target ---")
print(f"Total Fault Path Impedance (Z_total): {Z_total_pu.real:.3f} + j{Z_total_pu.imag:.3f} p.u.")
print(f"Power required by faulted line (P_W + jQ_W): {P_W_pu:.4f} + j{Q_W_pu:.4f} p.u.")
print(f"Real Power from Generator (P_gen, with losses): {P_gen_pu:.4f} p.u.")
print("\n--- Optimization based on Power Factor Constraint ---")
print(f"Max generator reactive power at {min_pf} PF (Q_gen_max): {Q_gen_max_pu:.4f} p.u.")
print("\n--- Final Equation for Optimal Reactive Power ---")
print(f"Q_opt (p.u.) = Q_W (p.u.) - Q_gen_max (p.u.)")
print(f"Q_opt (p.u.) = {Q_W_pu:.4f} - {Q_gen_max_pu:.4f}")
print(f"Q_opt (p.u.) = {Q_opt_pu:.4f} p.u.")
print("\n--- Final Answer ---")
print(f"The optimal reactive power injection required is: {Q_opt_MVAR:.4f} MVAR")

# Verification
if Q_opt_MVAR > Q_max_mvar:
    print(f"\nWarning: The required reactive power ({Q_opt_MVAR:.2f} MVAR) exceeds the compensator's maximum capacity of {Q_max_mvar} MVAR.")
else:
    print(f"\nThe calculated value is within the compensator's maximum capacity of {Q_max_mvar} MVAR.")

<<<0.2643>>>