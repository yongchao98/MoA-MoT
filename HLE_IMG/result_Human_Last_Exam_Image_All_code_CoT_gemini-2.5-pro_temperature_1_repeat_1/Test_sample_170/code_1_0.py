import cmath
import math

# --- Given Parameters ---
# Assuming impedances are in per-unit (p.u.) based on problem analysis
Z_S_pu = complex(0.02, 0.10)
Z_F_pu = complex(0.15, 0)
V_B_pre_comp_pu = 0.85  # Voltage at Bus B before compensation
V_B_target_pu = 1.0   # Target voltage at Bus B
S_base_MVA = 100
Q_max_MVAR = 50
harmonic_loss_increase = 0.04

# --- Step 1: Calculate the equivalent source voltage |V_A| ---
# In the pre-compensation state, the only load is the fault
I_F_pre_pu = V_B_pre_comp_pu / Z_F_pu
V_A_pu = V_B_pre_comp_pu + I_F_pre_pu * Z_S_pu
V_A_mag_pu = abs(V_A_pu)

# --- Step 2: Solve for the optimal reactive power Q_opt ---
# This requires solving the quadratic equation aQ^2 + bQ + c = 0 derived from the power flow equation.
# |V_A|^2 = |V_B_target + (I_F_comp + jQ) * Z_S|^2
I_F_comp_pu = V_B_target_pu / Z_F_pu

# Coefficients of the quadratic equation aQ^2 + bQ + c = 0
K1 = V_B_target_pu + I_F_comp_pu * Z_S_pu
K2 = complex(0, 1) * Z_S_pu

a = abs(K2)**2
b = 2 * (K1 * K2.conjugate()).real
c = abs(K1)**2 - V_A_mag_pu**2

# Solve the quadratic equation for Q
discriminant = b**2 - 4 * a * c

if discriminant < 0:
    print("Error: No real solution for Q_opt exists. The problem parameters may be inconsistent.")
    Q_opt_pu = 0
    is_solvable = False
else:
    # The minimum reactive power corresponds to the smaller valid root
    sol1 = (-b + math.sqrt(discriminant)) / (2 * a)
    sol2 = (-b - math.sqrt(discriminant)) / (2 * a)
    # The smaller positive root is the physically meaningful solution for minimum injection
    if sol2 > 0:
        Q_opt_pu = sol2
    else:
        Q_opt_pu = sol1
    is_solvable = True

# Convert to MVAR
Q_opt_MVAR = Q_opt_pu * S_base_MVA

# --- Step 3: Calculate system losses ---
# Line current in the compensated state
I_line_comp_pu = I_F_comp_pu + complex(0, Q_opt_pu)

# Real power loss in the transmission line Z_S
P_loss_pu = (abs(I_line_comp_pu)**2) * Z_S_pu.real

# Account for 4% harmonic losses
P_loss_total_pu = P_loss_pu * (1 + harmonic_loss_increase)

# Convert to MW
P_loss_total_MW = P_loss_total_pu * S_base_MVA

# --- Step 4: Print the final results ---
print("--- Problem Formulation and Solution ---")
print("\n1. System Parameters (in per-unit):")
print(f"   Source Impedance (Z_S): {Z_S_pu.real:.2f} + j{Z_S_pu.imag:.2f} p.u.")
print(f"   Fault Impedance (Z_F): {Z_F_pu.real:.2f} p.u.")
print(f"   Pre-fault Voltage at Bus B: {V_B_pre_comp_pu:.2f} p.u.")
print(f"   Base Power (S_base): {S_base_MVA} MVA")

print("\n2. Calculation of Equivalent Source Voltage |V_A|:")
print(f"   Source Voltage Magnitude |V_A| = {V_A_mag_pu:.3f} p.u.")

print("\n3. Optimization for Reactive Power Injection (Q_opt):")
print("   Objective: Minimize Q_comp, subject to V_B(Q_comp) = 1.0 p.u.")
print("   The final equation to solve for Q_comp is:")
print(f"   {a:.4f} * Q_comp^2 + ({b:.4f}) * Q_comp + ({c:.4f}) = 0")

if is_solvable:
    print(f"\n   Solving for Q_comp yields the minimum required reactive power.")
    print(f"   Optimal Reactive Power (Q_opt) = {Q_opt_pu:.3f} p.u.")
    print(f"   In physical units, Q_opt = {Q_opt_MVAR:.1f} MVAR")
    if Q_opt_MVAR > Q_max_MVAR:
        print(f"   Note: Required power ({Q_opt_MVAR:.1f} MVAR) exceeds STATCOM capacity ({Q_max_MVAR} MVAR). Result is based on problem's restoration requirement.")

print("\n4. Calculation of System Real Power Losses:")
print(f"   Total line current (I_line) = {I_line_comp_pu:.3f} p.u.")
print(f"   Line losses (P_loss) = |{abs(I_line_comp_pu):.3f}|^2 * {Z_S_pu.real:.2f} = {P_loss_pu:.3f} p.u.")
print(f"   Total System Losses with 4% harmonic effects = {P_loss_total_pu:.3f} p.u.")
print(f"   In physical units, Total System Losses = {P_loss_total_MW:.1f} MW")

# Final Answer
# print(f"<<<{Q_opt_MVAR:.1f}>>>")