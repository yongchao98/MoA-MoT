import math

def solve_hvac_optimization():
    """
    Calculates the optimal reactive power injection and system losses for the given HVAC system.
    """
    # --- Step 0: Define Givens and Base Values ---
    V_B_nom = 220e3  # V (Nominal voltage at Bus B)
    S_base = 100e6   # VA (Base power)
    V_B_fault_ratio = 0.85 # Voltage at Bus B during fault is 85% of nominal

    # Impedances (assuming values are given in per-unit, Ω is a typo)
    Z_S_pu = 0.02 + 1j * 0.10 # pu (Thevenin impedance Z_th = Z_S)
    Z_F_pu = 0.15 + 0j # pu (Fault impedance)
    R_th_pu = Z_S_pu.real
    X_th_pu = Z_S_pu.imag

    # STATCOM and System Constraints
    PF_min = 0.98
    harmonic_loss_increase = 0.04

    # --- Step 1: Characterize the system to match the fault condition ---
    # We find the Thevenin voltage (V_th) of the grid that results in |V_B| = 0.85 pu.
    V_B1_pu = V_B_fault_ratio
    P_L1_pu = V_B1_pu**2 / Z_F_pu.real
    Q_L1_pu = 0.0

    # Solve for |V_th|^2 using the full voltage equation:
    # |V_th|^2 = (|V_B| + (P*R + Q*X)/|V_B|)^2 + ((P*X - Q*R)/|V_B|)^2
    term1_Vth_sq = (V_B1_pu + (P_L1_pu * R_th_pu + Q_L1_pu * X_th_pu) / V_B1_pu)**2
    term2_Vth_sq = ((P_L1_pu * X_th_pu - Q_L1_pu * R_th_pu) / V_B1_pu)**2
    V_th_pu_sq = term1_Vth_sq + term2_Vth_sq
    V_th_pu = math.sqrt(V_th_pu_sq)

    # --- Step 2: Solve the Optimization Problem ---
    # Restore voltage to |V_B| = 1.0 pu while satisfying Bus B PF > 0.98.
    V_B2_pu = 1.0
    P_L2_pu = V_B2_pu**2 / Z_F_pu.real # Fault power at restored voltage

    # The minimum Q_comp occurs at the PF boundary: Q_comp = tan(acos(PF)) * (P_L2 - P_comp)
    tan_phi_max = math.tan(math.acos(PF_min))

    # This leads to a quadratic equation for Q_comp: a*Qc^2 + b*Qc + c = 0
    term_A = R_th_pu / tan_phi_max - X_th_pu
    term_B = X_th_pu / tan_phi_max + R_th_pu
    
    a = (term_A**2 + term_B**2) / V_B2_pu**2
    b = 2 * term_A
    c = V_B2_pu**2 - V_th_pu_sq

    # Solve the quadratic equation
    discriminant = b**2 - 4 * a * c
    # The valid solution for Q_comp is the smaller positive root
    sol1 = (-b + math.sqrt(discriminant)) / (2 * a)
    sol2 = (-b - math.sqrt(discriminant)) / (2 * a)
    Q_opt_pu = min(s for s in [sol1, sol2] if s > 0)
    Q_opt_MVAR = Q_opt_pu * S_base / 1e6

    # --- Step 3: Calculate System Losses ---
    # Find net power flow in the compensated state
    P_comp_pu = P_L2_pu - Q_opt_pu / tan_phi_max
    P_net_pu = P_L2_pu - P_comp_pu
    Q_net_pu = -Q_opt_pu
    S_net_pu = P_net_pu + 1j * Q_net_pu

    # Calculate losses: P_loss = |I_net|^2 * R_th * (1 + harmonic_factor)
    I_net_pu = abs(S_net_pu) / V_B2_pu
    P_loss_base_pu = I_net_pu**2 * R_th_pu
    P_loss_total_pu = P_loss_base_pu * (1 + harmonic_loss_increase)
    P_loss_total_MW = P_loss_total_pu * S_base / 1e6

    # --- Final Output ---
    print("--- Analysis and Final Equations ---")
    print("Based on the problem statement, a consistent model requires a Thevenin Voltage of "
          f"V_th = {V_th_pu:.4f} pu.")
    print("\nThe optimization problem is solved by substituting the power factor constraint "
          f"(Q_comp = {tan_phi_max:.4f} * (P_net)) into the voltage equation.")
    print("This results in the following quadratic equation for the optimal reactive power (Q_opt):")
    print(f"({a:.4f}) * Q_opt^2 + ({b:.4f}) * Q_opt + ({c:.4f}) = 0")
    print("\n--- Final Answer ---")
    print(f"The minimum reactive power injection required is Q_opt = {Q_opt_MVAR:.2f} MVAR.")
    print(f"The system's real power losses under these conditions are {P_loss_total_MW:.2f} MW.")

solve_hvac_optimization()
<<<The minimum reactive power injection required is 97.91 MVAR.The system's real power losses under these conditions are 50.34 MW.>>>