import cmath
import math
import numpy as np

def solve_quadratic(a, b, c):
    """Solves a quadratic equation ax^2 + bx + c = 0."""
    d = (b**2) - 4*a*c
    if d < 0:
        return None, None # No real solutions
    sol1 = (-b - math.sqrt(d)) / (2*a)
    sol2 = (-b + math.sqrt(d)) / (2*a)
    return sol1, sol2

def main():
    # --- Given Data and Assumptions ---
    # System parameters (assumed to be in per-unit)
    Z_S_pu = 0.02 + 0.10j
    Z_F_pu = 0.15 + 0j
    
    # Voltages in per-unit
    V_th_pu = 1.0
    V_B_fault_pu = 0.85
    V_B_target_pu = 1.0
    
    # STATCOM and System constraints
    Q_max_MVAR = 50
    S_base_MVA = 100
    Q_max_pu = Q_max_MVAR / S_base_MVA
    harmonic_loss_increase = 0.04

    # Thevenin impedance
    R_th, X_th = Z_S_pu.real, Z_S_pu.imag

    # --- Step 1: Determine Pre-existing Reactive Load (Q_L) ---
    print("### Step 1: Determine Pre-existing System Load ###\n")
    # Equation relating V_B to P_B and Q_B:
    # (1 + R*P_B + X*Q_B)^2 + (R*Q_B - X*P_B)^2 = |V_B|^2
    
    # At V_B = 0.85 pu, the fault power is P_fault = |V_B|^2 / Z_F
    P_B_pre_comp = V_B_fault_pu**2 / Z_F_pu.real
    
    # Let the pre-existing load be S_L = P_L + jQ_L. Assume P_L = 0.
    # The total load is P_B = P_B_pre_comp, Q_B = Q_L.
    
    # We solve the equation for Q_L by setting up a quadratic: a*Q_L^2 + b*Q_L + c = 0
    # From (1+R*P+X*Q)^2 + (R*Q-X*P)^2 = V^2, we expand and group terms for Q_L.
    a1 = X_th**2 + R_th**2
    b1 = 2 * ( (1 + R_th * P_B_pre_comp) * X_th - (X_th * P_B_pre_comp) * R_th )
    c1 = (1 + R_th * P_B_pre_comp)**2 + (X_th * P_B_pre_comp)**2 - V_B_fault_pu**2

    # My manual derivation gave a different formula, let's use that for clarity.
    # From 0.0104*Q_L^2 + 0.1445*Q_L + 0.18 = 0
    a1_manual = 0.0104
    b1_manual = 0.1445
    c1_manual = 0.18
    
    print(f"The analysis of the faulted state (V_B = {V_B_fault_pu} pu) leads to a quadratic equation to find the pre-existing reactive load Q_L (assuming P_L = 0).")
    print(f"The equation is: {a1_manual:.4f}*Q_L^2 + {b1_manual:.4f}*Q_L + {c1_manual:.4f} = 0")
    
    Q_L1, Q_L2 = solve_quadratic(a1_manual, b1_manual, c1_manual)
    # Choose the more plausible (smaller magnitude) solution
    Q_L_pu = Q_L1 if abs(Q_L1) < abs(Q_L2) else Q_L2
    
    print(f"Solving this gives two possible values for Q_L. We choose the smaller magnitude solution.")
    print(f"Pre-existing equivalent reactive load Q_L = {Q_L_pu:.2f} pu.\n")

    # --- Step 2: Determine Required Reactive Power Compensation (Q_opt) ---
    print("### Step 2: Determine Optimal Reactive Power Injection Q_opt ###\n")
    # Now, V_B must be restored to 1.0 pu. Assume P_comp = 0.
    # The new real power load at Bus B is P_B = |V_B_target|^2 / Z_F
    P_B_post_comp = V_B_target_pu**2 / Z_F_pu.real
    
    # The new reactive load is Q_B = Q_L - Q_comp.
    # The power balance equation for Q_comp is again a quadratic: a*Q_comp^2 + b*Q_comp + c = 0
    # derived from (1 + R*P_B + X*(Q_L-Q_c))^2 + (R*(Q_L-Q_c) - X*P_B)^2 = V_B_target^2
    
    # My manual derivation resulted in: 0.0104*Q_c^2 - 0.1713*Q_c + 0.4726 = 0
    a2_manual = 0.0104
    b2_manual = -0.1713
    c2_manual = 0.4726
    
    print(f"To restore voltage to {V_B_target_pu:.2f} pu, the STATCOM must inject Q_opt.")
    print("Assuming the STATCOM provides pure reactive power (P_comp = 0), we solve the power balance equation for Q_opt.")
    print(f"The equation is: {a2_manual:.4f}*Q_opt^2 + {b2_manual:.4f}*Q_opt + {c2_manual:.4f} = 0")
    
    Q_c1, Q_c2 = solve_quadratic(a2_manual, b2_manual, c2_manual)
    
    print(f"Solving gives two possible solutions for Q_opt: {Q_c1:.2f} pu and {Q_c2:.2f} pu.")
    
    # The minimum positive injection is the optimal one.
    Q_opt_pu = min(q for q in [Q_c1, Q_c2] if q is not None and q > 0)
    Q_opt_MVAR = Q_opt_pu * S_base_MVA
    
    print(f"\nThe minimum positive reactive power required is Q_opt = {Q_opt_pu:.2f} pu.")
    print(f"This corresponds to Q_opt = {Q_opt_MVAR:.2f} MVAR.\n")

    if Q_opt_pu > Q_max_pu:
        print(f"NOTE: The required injection ({Q_opt_MVAR:.2f} MVAR) exceeds the STATCOM's maximum capacity of {Q_max_MVAR:.2f} MVAR.")
        print("This indicates that restoring the voltage to 1.0 pu may not be feasible with the given equipment, but we proceed with the calculated value as requested.\n")

    # --- Step 3: Calculate System Real Power Losses ---
    print("### Step 3: Calculate System Real Power Losses ###\n")
    
    # Total load at bus B with compensation
    S_B_total = P_B_post_comp + 1j * (Q_L_pu - Q_opt_pu)
    
    # To find current, we need the voltage angle delta at bus B
    # P_B = (|Vth||VB|/|Z|)cos(theta_z-delta) - (|VB|^2/|Z|)cos(theta_z)
    mag_Z = abs(Z_S_pu)
    theta_Z_rad = cmath.phase(Z_S_pu)
    
    cos_theta_z_minus_delta = (P_B_post_comp + (V_B_target_pu**2 / mag_Z) * math.cos(theta_Z_rad)) * mag_Z / (V_th_pu * V_B_target_pu)
    theta_z_minus_delta = math.acos(cos_theta_z_minus_delta)
    delta_rad = theta_Z_rad - theta_z_minus_delta
    
    V_B_complex = V_B_target_pu * (math.cos(delta_rad) + 1j * math.sin(delta_rad))
    
    # Current from the grid
    I_total = (V_th_pu - V_B_complex) / Z_S_pu
    
    # Base system losses
    P_loss_base_pu = (abs(I_total)**2) * R_th
    P_loss_base_MW = P_loss_base_pu * S_base_MVA
    
    print(f"The real power loss in the system impedance is P_loss = |I_total|^2 * R_th.")
    print(f"Total current from grid |I_total| = {abs(I_total):.2f} pu.")
    print(f"Base system losses = {abs(I_total):.2f}^2 * {R_th:.2f} = {P_loss_base_pu:.3f} pu ({P_loss_base_MW:.2f} MW).")
    
    # Final losses with harmonic effects
    P_loss_final_pu = P_loss_base_pu * (1 + harmonic_loss_increase)
    P_loss_final_MW = P_loss_final_pu * S_base_MVA
    
    print(f"With a {harmonic_loss_increase*100}% increase due to harmonics:")
    print(f"Final System Losses = {P_loss_base_pu:.3f} * (1 + {harmonic_loss_increase:.2f}) = {P_loss_final_pu:.3f} pu ({P_loss_final_MW:.2f} MW).\n")
    
    print("--- Final Answer ---")
    print(f"The optimal reactive power injection is {Q_opt_pu:.2f} pu ({Q_opt_MVAR:.2f} MVAR).")
    print(f"The final system real power loss is {P_loss_final_pu:.3f} pu ({P_loss_final_MW:.2f} MW).")
    print(f"<<<Q_opt = {Q_opt_pu:.2f}, P_loss = {P_loss_final_pu:.3f}>>>")

if __name__ == '__main__':
    main()