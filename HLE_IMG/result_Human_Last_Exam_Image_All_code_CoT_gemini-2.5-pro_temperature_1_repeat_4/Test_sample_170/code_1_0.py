import cmath
import numpy as np

# --- Step 1: Define System Parameters in Per-Unit ---
S_base = 100.0  # MVA
V_base_B = 220.0 # kV
V_th_pu = 1.0 + 0j  # Thevenin voltage (external grid) in p.u.
Zs_pu = 0.02 + 0.1j # System impedance in p.u.
Vb_fault_pu = 0.85  # Voltage at Bus B during fault in p.u.
Vb_target_pu = 1.0  # Target voltage at Bus B after compensation in p.u.
harmonic_loss_factor = 1.04 # 4% increase in losses

# --- Step 2: Characterize the Fault Load ---
# We assume the fault is equivalent to a purely resistive shunt load R_fault.
# We find the value of R_fault that causes the voltage to drop to 0.85 p.u.
# The equation is |Vb_fault| = |V_th / (1 + Zs / R_fault)|
# Rearranging gives: |1 + Zs/R_fault|^2 = (1 / |Vb_fault|)^2
# Let x = 1/R_fault (the fault conductance G_fault)
# |(1 + Zs.real*x) + j(Zs.imag*x)|^2 = (1/Vb_fault_pu)**2
# (1 + Zs.real*x)^2 + (Zs.imag*x)^2 = (1/Vb_fault_pu)**2
# (1 + 2*Zs.real*x + Zs.real**2*x**2) + (Zs.imag**2*x**2) = (1/Vb_fault_pu)**2
# (Zs.real**2 + Zs.imag**2)*x**2 + (2*Zs.real)*x + (1 - (1/Vb_fault_pu)**2) = 0
# This is a quadratic equation for x = G_fault.

a = abs(Zs_pu)**2
b = 2 * Zs_pu.real
c = 1 - (1 / Vb_fault_pu)**2

# Solve the quadratic equation for x (G_fault)
discriminant = b**2 - 4 * a * c
if discriminant < 0:
    print("Error: No real solution for fault conductance exists.")
else:
    # We take the positive root for conductance
    G_fault_pu = (-b + np.sqrt(discriminant)) / (2 * a)
    R_fault_pu = 1 / G_fault_pu

    # --- Step 3: Solve for Optimal Reactive Power (Q_opt) ---
    # The STATCOM adds a shunt susceptance B_comp. Total admittance Y_total = G_fault + j*B_comp
    # We want |Vb_target| = |V_th / (1 + Zs * Y_total)| = 1.0
    # This means |1 + Zs * (G_fault + j*Q_comp)|^2 = 1.0 (since Q_comp = B_comp at 1 p.u. voltage)
    # Let Q = Q_comp
    # |1 + (Zs.real + j*Zs.imag)*(G_fault + j*Q)|^2 = 1.0
    # |(1 + Zs.real*G_fault - Zs.imag*Q) + j*(Zs.imag*G_fault + Zs.real*Q)|^2 = 1.0
    # (1 + Zs.real*G_fault - Zs.imag*Q)^2 + (Zs.imag*G_fault + Zs.real*Q)^2 = 1.0
    # This is a quadratic equation for Q.

    # Coefficients for aQ^2 + bQ + c = 0
    qa = abs(Zs_pu)**2
    qb = -2 * Zs_pu.imag * (1 + Zs_pu.real * G_fault_pu) + 2 * Zs_pu.real * (Zs_pu.imag * G_fault_pu)
    qb_simplified = -2 * Zs_pu.imag
    qc = (1 + Zs_pu.real * G_fault_pu)**2 + (Zs_pu.imag * G_fault_pu)**2 - 1.0

    # Solve the quadratic equation for Q_comp
    q_solutions = np.roots([qa, qb_simplified, qc])

    # The optimal Q is the minimum positive real solution
    positive_solutions = [s.real for s in q_solutions if s.real > 0]
    Q_opt_pu = min(positive_solutions)
    Q_opt_mvar = Q_opt_pu * S_base

    # --- Step 4: Calculate System Losses ---
    # Final admittance at Bus B
    Y_final_pu = G_fault_pu + 1j * Q_opt_pu
    
    # Final complex voltage at Bus B
    Vb_final_pu = V_th_pu / (1 + Zs_pu * Y_final_pu)
    
    # Line current from grid to Bus B
    I_line_pu = (V_th_pu - Vb_final_pu) / Zs_pu
    
    # System losses (real power loss in Zs)
    P_loss_pu = (abs(I_line_pu)**2) * Zs_pu.real
    
    # Add 4% for harmonic losses
    P_loss_final_pu = P_loss_pu * harmonic_loss_factor
    P_loss_final_mw = P_loss_final_pu * S_base

    # --- Step 5: Output the Results ---
    print("Optimization Results:")
    print(f"The minimum reactive power injection required (Q_opt) is {Q_opt_mvar:.2f} MVAR.")
    print(f"The final system real power loss is {P_loss_final_mw:.2f} MW.")
    print("\nCalculation Details:")
    print(f"Calculated fault conductance G_fault = {G_fault_pu:.4f} p.u.")
    print(f"Equation for Q_comp: ({qa:.4f})Q^2 + ({qb_simplified:.4f})Q + ({qc:.4f}) = 0")
    print(f"Solutions for Q_comp in p.u.: {q_solutions[0]:.4f}, {q_solutions[1]:.4f}")
    print(f"Chosen Q_opt = {Q_opt_pu:.4f} p.u.")
    print(f"Final line current = {abs(I_line_pu):.4f} p.u.")
    print(f"Base system losses = {P_loss_pu * S_base:.2f} MW")
    print(f"Final losses with 4% increase = {P_loss_final_mw:.2f} MW")
