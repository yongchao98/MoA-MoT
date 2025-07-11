import numpy as np
import cmath
from scipy.optimize import fsolve

def solve_hvac_optimization():
    """
    Solves the non-linear optimization problem for the HVAC transmission system.
    """
    # Step 1: System Parameters (p.u. and constants)
    S_base = 100.0  # MVA
    Z_S = 0.02 + 0.1j  # System impedance in p.u.
    Z_F = 0.15 + 0j    # Fault impedance in p.u.
    Q_max_statcom = 50.0  # MVAR
    Q_max_pu = Q_max_statcom / S_base  # 0.5 p.u.
    harmonic_loss_increase = 0.04

    # Admittances in p.u.
    Y_S = 1 / Z_S
    Y_F = 1 / Z_F
    G_S, B_S = Y_S.real, Y_S.imag
    G_F = Y_F.real

    # Step 2: Determine the optimal reactive power injection (Q_opt)
    # The optimization goal is to find the minimum Q_comp. The feasible range is
    # determined by the STATCOM's physical capacity (|Q_comp| <= Q_max_pu) and the
    # power factor constraint. The most negative Q_comp is limited by the capacity.
    Q_opt_pu = -Q_max_pu

    # Step 3: Find the system operating point (delta) for Q_opt
    # We need to solve the equation Q_comp(delta) = Q_opt_pu for the angle delta.
    # The formula for Q_comp is derived from nodal analysis:
    # Q_comp = (cos(d) * B_S + sin(d) * G_S) * |V_B|*|V_A| - |V_B|^2 * (B_S)
    # With V_A=1, V_B=1, this simplifies. Let's use the verified form from the thought process:
    # Q_comp(delta) = 9.615*cos(delta) + 1.923*sin(delta) - 9.615
    
    target_q = 9.615 + Q_opt_pu # The constant part moved to the right side

    def find_delta(d):
        """Function to find the root of, which is the angle delta."""
        return 9.615 * np.cos(d) + 1.923 * np.sin(d) - target_q

    # Solve for delta. An initial guess of -0.18 rad is used based on analysis.
    delta = fsolve(find_delta, -0.180)[0]

    # Step 4: Calculate the corresponding P_comp and system losses
    V_A = 1.0 + 0j
    V_B = cmath.rect(1.0, delta)  # Bus B voltage at 1.0 p.u. with angle delta

    # Calculate P_comp from its derived formula:
    # P_comp(delta) = 1.923*cos(delta) - 9.615*sin(delta) - 8.59
    P_comp_pu = 1.923 * np.cos(delta) - 9.615 * np.sin(delta) - (G_S + G_F)

    # Calculate individual power losses
    # Loss in transmission line Z_S
    I_S = (V_A - V_B) / Z_S
    P_loss_S_pu = (abs(I_S)**2) * Z_S.real

    # Loss in fault impedance Z_F
    I_F = V_B / Z_F
    P_loss_F_pu = (abs(I_F)**2) * Z_F.real

    # Loss in STATCOM (real power it consumes)
    P_loss_STATCOM_pu = -P_comp_pu

    # Sum losses to get fundamental loss
    P_loss_fundamental_pu = P_loss_S_pu + P_loss_F_pu + P_loss_STATCOM_pu

    # Add harmonic losses
    P_loss_total_pu = P_loss_fundamental_pu * (1 + harmonic_loss_increase)

    # Convert final p.u. values back to MW/MVAR
    Q_opt_MVAR = Q_opt_pu * S_base
    P_loss_total_MW = P_loss_total_pu * S_base
    P_loss_S_MW = P_loss_S_pu * S_base
    P_loss_F_MW = P_loss_F_pu * S_base
    P_loss_STATCOM_MW = P_loss_STATCOM_pu * S_base
    P_loss_fundamental_MW = P_loss_fundamental_pu * S_base

    # Step 5: Print the results
    print("--- Optimization Results ---")
    print(f"Optimal reactive power injection (Q_opt): {Q_opt_pu:.3f} p.u. = {Q_opt_MVAR:.1f} MVAR")
    print("\n--- System Power Losses ---")
    print(f"Total system real power loss: {P_loss_total_pu:.3f} p.u. = {P_loss_total_MW:.1f} MW")
    print("\n--- Loss Calculation Breakdown ---")
    print(f"Fundamental Loss = P_loss_line + P_loss_fault + P_loss_STATCOM")
    print(f"Fundamental Loss = {P_loss_S_MW:.1f} MW + {P_loss_F_MW:.1f} MW + {P_loss_STATCOM_MW:.1f} MW = {P_loss_fundamental_MW:.1f} MW")
    print(f"Total Loss = Fundamental Loss * (1 + Harmonic Increase)")
    print(f"Total Loss = {P_loss_fundamental_MW:.1f} MW * {1 + harmonic_loss_increase:.2f} = {P_loss_total_MW:.1f} MW")

solve_hvac_optimization()