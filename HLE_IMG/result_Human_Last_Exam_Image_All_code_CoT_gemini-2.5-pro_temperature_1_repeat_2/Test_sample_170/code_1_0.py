import cmath
import math

# --- Given Parameters ---
# Assuming these impedance values are in per-unit (p.u.) on the 100 MVA base.
Z_S_pu = 0.02 + 0.10j  # System impedance in p.u.
R_S_pu = Z_S_pu.real
X_S_pu = Z_S_pu.imag

Z_F_pu = 0.15 + 0.0j   # Fault impedance in p.u. (assumed resistive)

V_B_fault_pu = 0.85  # Voltage at Bus B during fault in p.u.
V_B_target_pu = 1.0  # Target voltage at Bus B in p.u.

S_base_MVA = 100.0   # Base power in MVA
Q_max_MVAR = 50.0    # Maximum STATCOM reactive power in MVAR
Q_max_pu = Q_max_MVAR / S_base_MVA

harmonic_loss_factor = 1.04 # 4% increase in losses

# --- Step 1: Calculate the Thevenin Voltage of the Grid ---
# Using the voltage divider formula for the uncompensated fault condition.
# V_B = V_th * (Z_F / (Z_S + Z_F)) => |V_th| = |V_B| * |(Z_S + Z_F) / Z_F|
# This is equivalent to |V_th| = |V_B| * |1 + Z_S / Z_F|

term = 1 + Z_S_pu / Z_F_pu
V_th_mag_pu = V_B_fault_pu * abs(term)

print("--- Step 1: System Analysis ---")
print(f"Ratio Z_S/Z_F = {Z_S_pu/Z_F_pu:.4f}")
print(f"Calculated Thevenin Voltage Magnitude |V_th| = {V_th_mag_pu:.4f} p.u.")
print("-" * 30)

# --- Step 2: Calculate the Required Reactive Power (Q_opt) ---
# Using the simplified real-valued power-voltage approximation to find Q_comp
# V_B^2 * (1 + R_S/Z_F) = V_B*|V_th| + Q_comp*X_S
# We solve for Q_comp when V_B is restored to 1.0 p.u.
# 1.0^2 * (1 + R_S/Z_F) = 1.0*|V_th| + Q_opt*X_S

lhs = V_B_target_pu**2 * (1 + R_S_pu / Z_F_pu.real)
rhs_V_part = V_B_target_pu * V_th_mag_pu

# Q_opt * X_S = lhs - rhs_V_part
Q_opt_pu = (lhs - rhs_V_part) / X_S_pu
Q_opt_MVAR = Q_opt_pu * S_base_MVA

print("--- Step 2: Optimal Reactive Power Calculation ---")
print(f"Objective: Restore voltage at Bus B to {V_B_target_pu:.2f} p.u.")
print(f"Required optimal reactive power Q_opt = {Q_opt_pu:.4f} p.u.")
print(f"This is equivalent to Q_opt = {Q_opt_MVAR:.2f} MVAR")
# Check against the STATCOM capacity
if Q_opt_MVAR > Q_max_MVAR:
    print(f"Warning: Required Q_opt ({Q_opt_MVAR:.2f} MVAR) exceeds Q_max ({Q_max_MVAR:.2f} MVAR).")
else:
    print(f"The required Q_opt is within the STATCOM's capacity of {Q_max_MVAR:.2f} MVAR.")
print("-" * 30)

# --- Step 3: Calculate System Real Power Losses ---
# Loss is primarily line loss P_loss = |I_AB|^2 * R_S
# Calculate current I_AB in the final compensated state (V_B = 1.0 p.u.)
# I_AB = I_fault + I_statcom = (V_B / Z_F) + (S_statcom* / V_B*)
# S_statcom = -jQ_opt (injection), S_statcom* = +jQ_opt
# Let's assume V_B is at angle 0 for simplicity in this calculation (V_B=1.0+0j)
V_B_final_pu = 1.0 + 0.0j

I_fault_pu = V_B_final_pu / Z_F_pu
I_statcom_pu = (1j * Q_opt_pu) / V_B_final_pu.conjugate()
I_AB_pu = I_fault_pu + I_statcom_pu

I_AB_mag_sq_pu = abs(I_AB_pu)**2

# Calculate base losses (line losses)
P_loss_base_pu = I_AB_mag_sq_pu * R_S_pu

# Apply the harmonic loss factor
P_loss_final_pu = P_loss_base_pu * harmonic_loss_factor
P_loss_final_MW = P_loss_final_pu * S_base_MVA

print("--- Step 3: System Loss Calculation ---")
print(f"Current from fault I_F = {I_fault_pu:.4f} p.u.")
print(f"Current from STATCOM I_STATCOM = {I_statcom_pu:.4f} p.u.")
print(f"Total line current I_AB = {I_AB_pu:.4f} p.u.")
print(f"Base line losses = |{I_AB_pu:.4f}|^2 * {R_S_pu:.4f} = {P_loss_base_pu:.4f} p.u.")
print(f"Final losses with 4% harmonic increase = {P_loss_base_pu:.4f} * {harmonic_loss_factor} = {P_loss_final_pu:.4f} p.u.")
print(f"This is equivalent to Final System Losses = {P_loss_final_MW:.2f} MW")
print("-" * 30)

# --- Final Answer Summary ---
print("\n--- Final Results ---")
print(f"The optimal reactive power injection required is Q_opt = {Q_opt_MVAR:.2f} MVAR.")
print(f"The final system real power losses are {P_loss_final_MW:.2f} MW.")

# Final equation for Q_opt
print("\nFinal Equation for Q_opt (in MVAR):")
print(f"Q_opt = [ (({V_B_target_pu:.2f}^2 * (1 + {R_S_pu:.2f}/{Z_F_pu.real:.2f})) - ({V_B_target_pu:.2f} * {V_th_mag_pu:.4f})) / {X_S_pu:.2f} ] * {S_base_MVA:.0f}")
print(f"Q_opt = [ (({V_B_target_pu**2 * (1 + R_S_pu/Z_F_pu.real):.4f}) - ({V_B_target_pu*V_th_mag_pu:.4f})) / {X_S_pu:.2f} ] * {S_base_MVA:.0f}")
print(f"Q_opt = [ {lhs - rhs_V_part:.4f} / {X_S_pu:.2f} ] * {S_base_MVA:.0f} = {Q_opt_MVAR:.2f} MVAR")

# Final equation for Power Loss
print("\nFinal Equation for Power Loss (in MW):")
print(f"P_loss = |({V_B_final_pu.real:.2f}/{Z_F_pu.real:.2f}) + j*({Q_opt_pu:.4f})|^2 * {R_S_pu:.2f} * {harmonic_loss_factor} * {S_base_MVA:.0f}")
print(f"P_loss = |{I_AB_pu:.4f}|^2 * {R_S_pu:.2f} * {harmonic_loss_factor} * {S_base_MVA:.0f}")
print(f"P_loss = {I_AB_mag_sq_pu:.4f} * {R_S_pu:.2f} * {harmonic_loss_factor} * {S_base_MVA:.0f} = {P_loss_final_MW:.2f} MW")