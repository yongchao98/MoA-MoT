import cmath
import math

def solve_hvac_system():
    """
    Solves the non-linear optimization problem for the HVAC transmission system.
    This function calculates the minimum reactive power injection from a STATCOM
    to restore bus voltage after a fault and computes the resulting system losses.
    """

    # --- 1. System Parameters (in per-unit) ---
    S_base = 100.0  # MVA
    V_B_nom = 220.0 # kV
    V_S = 1.0 + 0j # pu, Source voltage at external grid (assumed infinite bus)
    Z_S = 0.02 + 0.10j # pu, System impedance
    Z_F = 0.15 + 0j # pu, Fault impedance
    Q_max = 50.0 / S_base # pu, STATCOM max reactive power
    V_B_pre = 0.85 # pu, Voltage at Bus B during fault, pre-compensation
    V_B_post = 1.0 # pu, Target voltage at Bus B, post-compensation
    PF_min = 0.98 # STATCOM minimum power factor
    harmonic_loss_factor = 1.04

    # --- 2. Thevenin Equivalent Calculation ---
    Z_th = (Z_S * Z_F) / (Z_S + Z_F)
    R_th = Z_th.real
    X_th = Z_th.imag
    
    # --- 3. & 4. & 5. Solve for Optimal Reactive Power (Q_opt) ---
    # From PF > 0.98, we get |P_c| / sqrt(P_c^2 + Q_c^2) > 0.98
    # which implies |P_c| > 4.92 * |Q_c|.
    # We use the boundary P_c = 4.92 * Q_c to minimize Q_c.
    pf_ratio = math.sqrt(1/PF_min**2 - 1)
    Pc_over_Qc = 1 / pf_ratio # P_c = (1/pf_ratio) * Q_c
    
    # Use linearized voltage sensitivity: delta(|V|^2) = 2 * (P_c*R_th + Q_c*X_th)
    delta_V_sq = V_B_post**2 - V_B_pre**2

    # Substitute P_c = Pc_over_Qc * Q_c into the sensitivity equation
    Q_opt_pu = delta_V_sq / (2 * (Pc_over_Qc * R_th + X_th))
    Q_opt_mvar = Q_opt_pu * S_base

    # Required real power from STATCOM
    P_c_pu = Pc_over_Qc * Q_opt_pu
    P_c_mw = P_c_pu * S_base
    
    print("--- Optimization Results ---")
    print(f"Thevenin Impedance Z_th = ({R_th:.4f} + j{X_th:.4f}) pu")
    print(f"Required Voltage Change delta(|V|^2) = {V_B_post**2} - {V_B_pre**2} = {delta_V_sq:.4f} pu")
    print("Optimization equation: delta(|V|^2) = 2 * (P_c * R_th + Q_c * X_th)")
    print(f"Using the power factor constraint, P_c = {Pc_over_Qc:.2f} * Q_c. The equation becomes:")
    print(f"{delta_V_sq:.4f} = 2 * Q_c * ({Pc_over_Qc:.2f} * {R_th:.4f} + {X_th:.4f})")
    print(f"Optimal Reactive Power Injection (Q_opt) = {Q_opt_mvar:.2f} MVAR")
    print(f"Required Real Power Injection (P_c) = {P_c_mw:.2f} MW")

    # --- 6. & 7. Calculate System Losses ---
    # Loss calculation requires the final voltage angle, which is unknown without the load profile.
    # We assume a reasonable small angle of -5 degrees for a loaded line exporting power.
    delta_2_deg = -5.0
    delta_2_rad = math.radians(delta_2_deg)
    V_B_final = cmath.rect(V_B_post, delta_2_rad)

    # Line Loss
    R_S = Z_S.real
    I_S = (V_S - V_B_final) / Z_S
    P_loss_line_pu = (abs(I_S)**2) * R_S

    # Fault Loss (assuming Z_F is resistive)
    P_loss_fault_pu = (abs(V_B_final)**2) / Z_F.real

    # Total fundamental loss
    P_loss_total_pu = P_loss_line_pu + P_loss_fault_pu
    
    # Total loss including harmonics
    P_loss_final_pu = P_loss_total_pu * harmonic_loss_factor
    P_loss_final_mw = P_loss_final_pu * S_base

    print("\n--- System Loss Calculation (with assumed final angle) ---")
    print(f"Assumed final voltage angle at Bus B: {delta_2_deg} degrees")
    print(f"Line Loss Equation: |(V_S - V_B) / Z_S|^2 * R_S")
    print(f"Line Loss = |({V_S.real:.2f}+j{V_S.imag:.2f}) - ({V_B_final.real:.2f}+j{V_B_final.imag:.2f}) / ({Z_S.real:.2f}+j{Z_S.imag:.2f})|^2 * {R_S:.2f} = {P_loss_line_pu * S_base:.2f} MW")
    print(f"Fault Loss Equation: |V_B|^2 / R_F")
    print(f"Fault Loss = |{abs(V_B_final):.2f}|^2 / {Z_F.real:.2f} = {P_loss_fault_pu * S_base:.2f} MW")
    print(f"Total Fundamental Power Loss = {P_loss_line_pu*S_base:.2f} + {P_loss_fault_pu*S_base:.2f} = {P_loss_total_pu * S_base:.2f} MW")
    print(f"Total System Loss with 4% harmonics = {P_loss_total_pu * S_base:.2f} MW * {harmonic_loss_factor} = {P_loss_final_mw:.2f} MW")
    
    print(f"\n--- Final Answer ---")
    print(f"The optimal reactive power injection is {Q_opt_mvar:.1f} MVAR.")
    print(f"The system's real power losses are {P_loss_final_mw:.1f} MW.")
    
    return Q_opt_mvar

final_answer = solve_hvac_system()
print(f'<<<{final_answer:.1f}>>>')