import math

def solve_hvac_optimization():
    """
    Solves the non-linear optimization problem for the HVAC transmission system.
    """

    # --- Step 1: Define constants and assumptions in per-unit (p.u.) ---
    S_base = 100.0  # MVA
    V_B_nom = 220.0  # kV
    V_A = 1.0  # p.u., Source voltage at Bus A
    
    # Impedances are assumed to be in p.u. based on the context
    R_S = 0.02  # p.u.
    X_S = 0.10  # p.u.
    
    # Voltages at Bus B
    V_B_fault = 0.85  # p.u.
    V_B_final = 1.0  # p.u.
    
    # STATCOM maximum capacity
    Q_max_mvar = 50.0 # MVAR
    Q_max_pu = Q_max_mvar / S_base # p.u.
    
    # Harmonic loss factor
    harmonic_loss_factor = 1.04

    # --- Step 2 & 3: Calculate the optimal reactive power Q_opt ---
    # The voltage-power relationship is: V_B^2 = V_A^2 - 2*(R*P + X*Q)
    #
    # Before compensation:
    # V_B_fault^2 = V_A^2 - 2 * (R_S*P_L + X_S*Q_L)  (Eq. 1)
    #
    # After compensation (restoring V_B to 1.0 p.u.):
    # V_B_final^2 = V_A^2 - 2 * (R_S*P_L + X_S*(Q_L - Q_opt))  (Eq. 2)
    #
    # From Eq. 2, since V_B_final = V_A = 1.0:
    # 1.0^2 = 1.0^2 - 2 * (R_S*P_L + X_S*(Q_L - Q_opt))
    # 0 = R_S*P_L + X_S*Q_L - X_S*Q_opt
    # R_S*P_L + X_S*Q_L = X_S*Q_opt
    #
    # Substitute this back into Eq. 1:
    # V_B_fault^2 = V_A^2 - 2 * (X_S * Q_opt)
    #
    # Rearranging for Q_opt:
    # 2 * X_S * Q_opt = V_A^2 - V_B_fault^2
    # Q_opt = (V_A^2 - V_B_fault^2) / (2 * X_S)

    Q_opt_pu = (V_A**2 - V_B_fault**2) / (2 * X_S)
    Q_opt_mvar = Q_opt_pu * S_base
    
    print("--- Calculation of Optimal Reactive Power (Q_opt) ---")
    print(f"Equation for Q_opt: (V_A^2 - V_B_fault^2) / (2 * X_S)")
    print(f"Q_opt = ({V_A**2:.2f} - {V_B_fault**2:.4f}) / (2 * {X_S})")
    print(f"Q_opt = {V_A**2 - V_B_fault**2:.4f} / {2 * X_S:.2f}")
    print(f"Optimal reactive power (Q_opt) = {Q_opt_pu:.4f} p.u.")
    print(f"Optimal reactive power (Q_opt) = {Q_opt_mvar:.2f} MVAR\n")

    # --- Step 4: Calculate system losses ---
    # The real power load (P_L) is not given. We must assume a value to calculate losses.
    # A common assumption is to use the base power for the load.
    P_L_pu = 1.0  # Assume load is 100 MW = 1.0 p.u.
    
    # From the derivation in Step 3 (Eq. 2), with V_B restored to 1.0:
    # 0 = R_S*P_L + X_S*(Q_L - Q_opt)
    # Let Q'_L = Q_L - Q_opt (the net reactive load after compensation)
    # 0 = R_S*P_L + X_S*Q'_L
    # Q'_L = -(R_S / X_S) * P_L
    Q_L_final_pu = -(R_S / X_S) * P_L_pu
    
    # Calculate the line current after compensation
    # |I_line|^2 = P_L^2 + Q'_L^2  (since V_B_final = 1.0 p.u.)
    I_line_sq_pu = P_L_pu**2 + Q_L_final_pu**2
    
    # Calculate line losses
    P_loss_line_pu = I_line_sq_pu * R_S
    
    # Calculate total system losses including harmonic effects
    P_loss_total_pu = P_loss_line_pu * harmonic_loss_factor
    P_loss_total_mw = P_loss_total_pu * S_base
    
    print("--- Calculation of System Real Power Losses ---")
    print(f"Assumption: Active Power Load (P_L) = {P_L_pu * S_base:.1f} MW ({P_L_pu:.2f} p.u.)")
    print(f"Net Reactive Power after compensation Q'_L = - (R_S / X_S) * P_L = -({R_S}/{X_S}) * {P_L_pu:.2f} = {Q_L_final_pu:.4f} p.u.")
    print(f"Line current squared |I'|^2 = P_L^2 + Q'_L^2 = {P_L_pu**2:.4f} + {Q_L_final_pu**2:.4f} = {I_line_sq_pu:.4f} p.u.^2")
    print(f"Line losses P_loss_line = |I'|^2 * R_S = {I_line_sq_pu:.4f} * {R_S} = {P_loss_line_pu:.4f} p.u.")
    print(f"Total losses P_loss_total = P_loss_line * (1 + 0.04) = {P_loss_line_pu:.4f} * {harmonic_loss_factor} = {P_loss_total_pu:.4f} p.u.")
    print(f"Total system real power losses = {P_loss_total_mw:.2f} MW")

# Execute the function
solve_hvac_optimization()