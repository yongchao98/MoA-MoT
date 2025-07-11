import numpy as np
import math

def solve_power_system_optimization():
    """
    Solves the power system optimization problem described.
    It calculates the required reactive power from a STATCOM to restore voltage after a fault
    and determines the resulting system power losses.
    """
    # --- Given System Parameters ---
    # Note: Impedances Z_S and Z_F are given in Ohms, but their values are typical for per-unit (p.u.).
    # Using them as Ohm values with the given base values results in extremely low p.u. impedances
    # and physically inconsistent results (e.g., a voltage sag to 0.85 p.u. is impossible).
    # We will proceed under the common assumption for such problems that 'Ω' is a typo for 'p.u.'.

    # Thevenin equivalent impedance of the grid as seen from the fault point/Bus B
    Z_th_pu = 0.02 + 0.10j
    R_th = Z_th_pu.real
    X_th = Z_th_pu.imag

    # Fault impedance (assumed to be purely resistive and in p.u.)
    Z_F_pu = 0.15

    # System voltages in per-unit (p.u.)
    V_th_pu = 1.0  # Thevenin equivalent voltage of the external grid
    V_B_fault_pu = 0.85  # Voltage at Bus B during the fault, before compensation
    V_B_target_pu = 1.0  # Target voltage at Bus B after compensation

    # System base power and harmonic loss factor
    S_base_MVA = 100.0
    harmonic_loss_increase_factor = 0.04

    # --- Step 1: Analyze pre-compensation state to find wind farm power P_W ---
    # In the state with the fault but without STATCOM compensation, we assume the wind farm
    # operates at unity power factor (Q_W = 0). The fault is resistive (Q_F = 0).
    # For power balance, the reactive power supplied from the grid to Bus B must be zero (Q_grid1 = 0).
    
    # The condition Q_grid1 = 0 leads to the equation:
    # V_B_fault * X_th = X_th*cos(delta_1) + R_th*sin(delta_1)
    # 0.85 * 0.10 = 0.10*cos(delta_1) + 0.02*sin(delta_1)
    # 0.085 = 0.10*cos(delta_1) + 0.02*sin(delta_1)
    
    # Solve for delta_1. We use an analytical solution method.
    # a*cos(x) + b*sin(x) = c  =>  sqrt(a^2+b^2)*cos(x-atan(b/a)) = c
    a, b, c = 0.10, 0.02, 0.085
    # The angle delta_1 is expected to be negative as power flows from the grid (angle 0) to Bus B.
    delta_1_rad = math.atan2(b, a) - math.acos(c / math.sqrt(a**2 + b**2))
    delta_1_deg = math.degrees(delta_1_rad)

    # Now calculate the real power from the grid (P_grid1)
    # Power flow from Bus B to the grid:
    # P_grid = (|V_B|^2 - |V_B||V_th|cos(d))R + (|V_B||V_th|sin(d))X / (R^2+X^2)
    # Q_grid = (|V_B|^2 - |V_B||V_th|cos(d))X - (|V_B||V_th|sin(d))R / (R^2+X^2)
    
    V_B1 = V_B_fault_pu
    d1 = delta_1_rad
    Z_th_sq_mag = R_th**2 + X_th**2
    
    P_grid1_num = (V_B1**2 - V_B1*V_th_pu*math.cos(d1))*R_th + (V_B1*V_th_pu*math.sin(d1))*X_th
    P_grid1_pu = P_grid1_num / Z_th_sq_mag
    
    # Calculate real power consumed by the fault (P_F1)
    P_F1_pu = V_B1**2 / Z_F_pu
    
    # Calculate wind farm power injection (P_W) from power balance at Bus B
    # P_W = P_F1 + P_grid1
    P_W_pu = P_F1_pu + P_grid1_pu

    # --- Step 2: Analyze post-compensation state to find required Q_opt ---
    # With STATCOM, the voltage is restored to V_B_target_pu = 1.0.
    # P_W is assumed constant.

    # Calculate real power consumed by the fault at target voltage (P_F2)
    V_B2 = V_B_target_pu
    P_F2_pu = V_B2**2 / Z_F_pu
    
    # Calculate real power from the grid (P_grid2) from power balance
    # P_W = P_F2 + P_grid2
    P_grid2_pu = P_W_pu - P_F2_pu
    
    # Solve for the new bus angle (delta_2) using the P_grid2 equation.
    # P_grid2 = ((V_B2^2 - V_B2*V_th*cos(d2))R + (V_B2*V_th*sin(d2))X) / (R^2+X^2)
    # A*sin(d2) - B*cos(d2) = C
    # Where A=V_B2*V_th*X, B=V_B2*V_th*R, C = P_grid2*(R^2+X^2) - V_B2^2*R
    A = V_B2 * V_th_pu * X_th
    B = V_B2 * V_th_pu * R_th
    C = P_grid2_pu * Z_th_sq_mag - V_B2**2 * R_th
    # Solving A*sin(x) - B*cos(x) = C is equivalent to sqrt(A^2+B^2)*sin(x-atan(B/A)) = C
    delta_2_rad = math.atan2(B, A) + math.asin(C / math.sqrt(A**2 + B**2))
    delta_2_deg = math.degrees(delta_2_rad)

    # Now calculate the reactive power from the grid (Q_grid2) with V_B2 and delta_2
    d2 = delta_2_rad
    Q_grid2_num = (V_B2**2 - V_B2*V_th_pu*math.cos(d2))*X_th - (V_B2*V_th_pu*math.sin(d2))*R_th
    Q_grid2_pu = Q_grid2_num / Z_th_sq_mag
    
    # The required STATCOM reactive power Q_opt is found from the reactive power balance:
    # Q_opt_injected = Q_consumed = Q_grid2 + Q_F2 - Q_W
    # Since Q_F2 = 0 (resistive fault) and Q_W = 0 (unity pf wind farm), Q_opt = Q_grid2
    Q_opt_pu = Q_grid2_pu
    Q_opt_MVAR = Q_opt_pu * S_base_MVA
    
    # --- Step 3: Calculate System Losses ---
    # System losses are the real power losses in the grid impedance Z_th.
    # P_loss = P_from_source - P_delivered_to_bus = P_th,B - P_B,th
    # We already calculated P_grid2, which is P_B,th. Let's calculate P_th,B.
    # P_th,B = |V_th|^2*G - |V_th||V_B2|*(G*cos(d_th-d_B) + B*sin(d_th-d_B))
    # where Y_th = G+jB = 1/Z_th
    Y_th = 1 / Z_th_pu
    G = Y_th.real
    B = Y_th.imag
    
    P_th_B_pu = V_th_pu**2*G - V_th_pu*V_B2*(G*math.cos(0-d2) + B*math.sin(0-d2))
    P_loss_pu = P_th_B_pu + P_grid2_pu # P_grid2 is power from B to th, so add it
    
    # Add the 4% increase due to harmonic distortion
    P_loss_total_pu = P_loss_pu * (1 + harmonic_loss_increase_factor)
    P_loss_total_MW = P_loss_total_pu * S_base_MVA

    # --- Print Results ---
    print("--- Problem Analysis and Results ---\n")
    print("This script solves for the required reactive power compensation and system losses")
    print("under the specified fault conditions.\n")
    
    print("Part 1: Required Reactive Power (Q_opt)\n")
    print("The optimal reactive power Q_opt is determined by solving the power flow equations")
    print(f"to achieve a target voltage of {V_B_target_pu:.2f} p.u. at Bus B.")
    print("\nThe equation for the required reactive power injection is Q_opt = Q_grid, where:")
    print("Q_grid = ((|V_B|^2 - |V_B|*|V_th|*cos(delta))*X_th - (|V_B|*|V_th|*sin(delta))*R_th) / (R_th^2 + X_th^2)\n")
    print("Using the calculated values:")
    print(f"  |V_B|    = {V_B2:.2f} p.u.")
    print(f"  |V_th|   = {V_th_pu:.2f} p.u.")
    print(f"  delta    = {delta_2_deg:.2f} degrees")
    print(f"  X_th     = {X_th:.2f} p.u.")
    print(f"  R_th     = {R_th:.2f} p.u.")
    print(f"  Q_opt    = (({V_B2:.2f}^2 - {V_B2:.2f}*{V_th_pu:.2f}*cos({delta_2_deg:.2f}))*{X_th:.2f} - ({V_B2:.2f}*{V_th_pu:.2f}*sin({delta_2_deg:.2f}))*{R_th:.2f}) / ({R_th:.2f}^2 + {X_th:.2f}^2)")
    print(f"  Q_opt    = {Q_opt_pu:.4f} p.u.\n")
    print(f"The minimum reactive power injection required from the MMCC STATCOM is:")
    print(f"Q_opt = {Q_opt_MVAR:.2f} MVAR")
    
    print("\n" + "-"*40 + "\n")
    
    print("Part 2: System Real Power Losses\n")
    print("System losses are calculated as the real power dissipated in the Thevenin impedance Z_th,")
    print("with an additional 4% to account for harmonic effects.")
    print("\nThe equation for the system losses is P_loss = P_loss_base * (1 + harmonic_factor), where P_loss_base can be found from the power mismatch across the line.\n")
    print("Calculated base loss P_loss_base:")
    print(f"  P_loss_base = {P_loss_pu:.4f} p.u.")
    print("Applying the harmonic loss factor:")
    print(f"  P_loss_total = P_loss_base * (1 + {harmonic_loss_increase_factor})")
    print(f"  P_loss_total = {P_loss_pu:.4f} * {1 + harmonic_loss_increase_factor}")
    print(f"  P_loss_total = {P_loss_total_pu:.4f} p.u.\n")
    print(f"The total system real power losses are:")
    print(f"P_loss = {P_loss_total_MW:.2f} MW")

# Execute the function
solve_power_system_optimization()