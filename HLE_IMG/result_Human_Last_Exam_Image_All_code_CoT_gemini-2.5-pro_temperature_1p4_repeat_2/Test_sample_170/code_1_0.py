import numpy as np
from scipy.optimize import fsolve

# Step 1: Define System Parameters
S_base = 100.0  # MVA
V_base_B = 220.0  # kV
V_b_fault_mag_pu = 0.85  # per unit
V_b_target_mag_pu = 1.0  # per unit

# Per-unit impedances
# Assuming the given impedance values are already in per unit on the system base.
Z_S_pu = 0.02 + 1j * 0.10
Z_F_pu = 0.15 + 1j * 0.0

# STATCOM maximum capacity
Q_max = 50.0 # MVAR
Q_max_pu = Q_max / S_base

# Harmonic loss factor
harmonic_loss_increase = 0.04

# Step 2: Model the System and State Assumptions
# We assume TR1 is ideal, so Bus A and Bus B are the same node in pu.
# The fault is modeled at Bus B. The STATCOM injects power at Bus B.

# Step 3: Calculate Thevenin Equivalent Impedance at Bus B
# Z_th is the parallel combination of Z_S and Z_F
Z_th_pu = (Z_S_pu * Z_F_pu) / (Z_S_pu + Z_F_pu)
R_th_pu = Z_th_pu.real
X_th_pu = Z_th_pu.imag

# Step 4: Formulate and Solve the Non-Linear Voltage Equation
# We need to solve for Q_opt in the equation:
# |V_target| - |V_fault_phasor + I_comp*Z_th| = 0
# V_target is 1.0 at angle 0. Let V_fault = 0.85 at angle delta.
# 1.0 - (0.85*cos(delta) + j*0.85*sin(delta)) = Q_opt * X_th - j*Q_opt*R_th
# From this we derive the non-linear equation for Q_opt to solve:
# 1.0 - 0.85 * sqrt(1 - (Q*R_th/0.85)^2) = Q * X_th

def equation_to_solve(Q_pu, R_th, X_th, V_fault, V_target):
    """
    Non-linear equation for reactive power compensation Q_pu.
    Derived from the exact phasor relationship.
    """
    # Argument of the square root, must be non-negative
    arg_sqrt = 1 - (Q_pu * R_th / V_fault)**2
    if arg_sqrt < 0:
        # Return a large number if Q is out of the valid domain
        return 1e6
    
    # Left side of the equation: V_target - V_fault_real_part
    lhs = V_target - V_fault * np.sqrt(arg_sqrt)
    # Right side of the equation: Q_pu * X_th
    rhs = Q_pu * X_th
    return lhs - rhs

# Solve for Q_opt_pu numerically
# An initial guess is made from the linear approximation: Q ~ (V_target-V_fault)/X_th
initial_guess_q = (V_b_target_mag_pu - V_b_fault_mag_pu) / X_th_pu
Q_opt_pu = fsolve(equation_to_solve, x0=initial_guess_q, args=(R_th_pu, X_th_pu, V_b_fault_mag_pu, V_b_target_mag_pu))[0]

# Convert optimal reactive power to MVAR
Q_opt_mvar = Q_opt_pu * S_base

# Step 5: Calculate System Power Losses
# After compensation, assume V_B is restored to 1.0 pu
V_B_final_pu = 1.0 + 0j

# Current into the fault
I_F_pu = V_B_final_pu / Z_F_pu

# Current injected by the STATCOM
# I_comp = (S_comp/V_B_final)*, S_comp = jQ_opt
I_comp_pu = (1j * Q_opt_pu / V_B_final_pu).conjugate()

# Current from the source (Grid) assuming no other load/generation
I_S_pu = I_F_pu - I_comp_pu

# Power loss in the source impedance Z_S
P_loss_S_pu = (abs(I_S_pu)**2) * Z_S_pu.real

# Power loss in the fault impedance Z_F
P_loss_F_pu = (abs(I_F_pu)**2) * Z_F_pu.real

# Total losses before considering harmonics
total_loss_pu = P_loss_S_pu + P_loss_F_pu

# Total losses including 4% increase due to harmonics
final_total_loss_pu = total_loss_pu * (1 + harmonic_loss_increase)

# Convert final losses to MW
final_total_loss_mw = final_total_loss_pu * S_base

# Step 6: Output the results
print("--- Optimization Results ---")
print(f"The calculated Thevenin impedance at Bus B is Z_th = {R_th_pu:.4f} + j{X_th_pu:.4f} pu.")
print("\n--- Required Reactive Power Compensation ---")
print(f"The equation solved for optimal reactive power (Q_opt) is:")
print(f"{V_b_target_mag_pu:.2f} - {V_b_fault_mag_pu:.2f} * sqrt(1 - (Q_opt * {R_th_pu:.4f} / {V_b_fault_mag_pu:.2f})^2) = Q_opt * {X_th_pu:.4f}")
print(f"Solving this gives Q_opt = {Q_opt_pu:.2f} pu.")
print(f"The minimum reactive power injection required from the STATCOM is Q_opt = {Q_opt_mvar:.2f} MVAR.")
if Q_opt_mvar > Q_max:
    print(f"(Note: This required value of {Q_opt_mvar:.2f} MVAR exceeds the STATCOM's maximum capacity of {Q_max:.1f} MVAR.)")

print("\n--- System Real Power Loss Calculation ---")
print("After compensation, with V_B restored to 1.0 pu:")
print(f"Final total system loss (including 4% harmonic effects) is calculated as:")
print(f"Loss = (1 + {harmonic_loss_increase}) * (|I_S|^2 * R_S + |I_F|^2 * R_F)")
print(f"Loss = (1 + {harmonic_loss_increase:.2f}) * (({abs(I_S_pu):.2f})^2 * {Z_S_pu.real:.2f} + ({abs(I_F_pu):.2f})^2 * {Z_F_pu.real:.2f})")
print(f"Loss = {1 + harmonic_loss_increase:.2f} * ({P_loss_S_pu:.2f} + {P_loss_F_pu:.2f}) pu = {final_total_loss_pu:.2f} pu.")
print(f"The total system real power loss is {final_total_loss_mw:.2f} MW.")
print("\n<<<Q_opt_mvar={:.2f},Loss_mw={:.2f}>>>".format(Q_opt_mvar, final_total_loss_mw))