import cmath
import math

def solve_hvac_system():
    # Step 1: Given data and Per-Unit conversion
    S_base = 100.0  # MVA
    V_B_nom = 220.0  # kV
    
    # Assume given impedances are in p.u.
    Z_S_pu = complex(0.02, 0.10)
    Z_F_pu = complex(0.15, 0)
    
    # Assume a typical transformer impedance for TR1
    Z_TR1_pu = complex(0, 0.20)
    
    # STATCOM constraints
    Q_max_mvar = 50.0
    Q_max_pu = Q_max_mvar / S_base

    # Voltages
    V_B_sag_pu = 0.85
    V_B_target_pu = 1.0
    
    # Harmonic loss factor
    harmonic_loss_factor = 1.04

    print("Step 1: System Parameters (in p.u.)")
    print(f"S_base = {S_base} MVA")
    print(f"Source Impedance Z_S = {Z_S_pu:.4f} p.u.")
    print(f"Fault Impedance Z_F = {Z_F_pu:.4f} p.u.")
    print(f"Transformer TR1 Impedance Z_TR1 = {Z_TR1_pu:.4f} p.u.")
    print(f"Voltage at Bus B during fault = {V_B_sag_pu:.2f} p.u.")
    print("-" * 30)

    # Step 2: Fault Modeling using Symmetrical Components
    Z1 = Z_S_pu
    Z2 = Z1  # Assuming Z2 = Z1
    Z0 = 3 * Z1 # Assuming Z0 = 3*Z1
    
    # Thevenin impedance at Bus A (fault location) for a SLG fault
    Z_th_A = (Z1 * (Z2 + Z0 + 3*Z_F_pu)) / (Z1 + Z2 + Z0 + 3*Z_F_pu)

    # Thevenin impedance at Bus B
    Z_th_B = Z_th_A + Z_TR1_pu
    R_th = Z_th_B.real
    X_th = Z_th_B.imag
    
    print("Step 2: Thevenin Equivalent during Fault")
    print(f"Thevenin Impedance at Bus B, Z_th_B = {Z_th_B.real:.4f} + j{Z_th_B.imag:.4f} p.u.")

    # Thevenin voltage at Bus B during fault
    V_th_B_mag = V_B_sag_pu
    print(f"Thevenin Voltage Magnitude at Bus B, |V_th_B| = {V_th_B_mag:.4f} p.u.")
    print("-" * 30)
    
    # Step 3: Solve for optimal reactive power Q_opt
    # We solve the equation set derived from V'_B * V'_B_conj = V_th_B * V'_B_conj + j*Q_c*Z_th_B
    # |V'_B|^2 = |V_th_B|*|V'_B|*cos(delta_th - delta_B) - Q_c*X_th
    # 0 = |V_th_B|*|V'_B|*sin(delta_th - delta_B) + Q_c*R_th
    # Let V_th_B be the reference (delta_th = 0), and |V'_B|=1.0
    # 1.0 = V_th_B_mag * cos(-delta_B) - Q_c * X_th
    # 0 = V_th_B_mag * sin(-delta_B) + Q_c * R_th => Q_c = V_th_B_mag * sin(delta_B) / R_th

    # Solve for delta_B
    # 1.0 = V_th_B_mag * cos(delta_B) - (V_th_B_mag * sin(delta_B) / R_th) * X_th
    # 1.0 / V_th_B_mag = cos(delta_B) - (X_th / R_th) * sin(delta_B)
    A = 1.0
    B = -X_th / R_th
    C = 1.0 / V_th_B_mag
    
    R_solver = math.sqrt(A**2 + B**2)
    alpha = math.atan2(B, A)
    
    # R_solver * cos(delta_B + alpha) = C => cos(delta_B + alpha) = C / R_solver
    cos_val = C / R_solver
    # Handling potential math domain error if system cannot be restored to 1.0 p.u.
    if abs(cos_val) > 1:
        print("Error: Voltage cannot be restored to 1.0 p.u. with the given system parameters.")
        return
        
    delta_B_plus_alpha = math.acos(cos_val)
    
    # Two possible solutions for delta_B, choose the smaller one
    delta_B1 = delta_B_plus_alpha - alpha
    delta_B2 = -delta_B_plus_alpha - alpha
    
    delta_B = min(abs(delta_B1), abs(delta_B2))

    # Calculate Q_opt in p.u.
    Q_opt_pu = (V_th_B_mag * math.sin(delta_B)) / R_th
    Q_opt_mvar = Q_opt_pu * S_base
    
    print("Step 3: Required Reactive Power Calculation")
    print(f"STATCOM voltage angle delta_B = {math.degrees(delta_B):.4f} degrees")
    print(f"Optimal Reactive Power Q_opt = {Q_opt_pu:.4f} p.u.")
    print(f"This is equivalent to Q_opt = {Q_opt_mvar:.4f} MVAR")

    if Q_opt_mvar > Q_max_mvar:
        print(f"Warning: Required Q_opt ({Q_opt_mvar:.2f} MVAR) exceeds STATCOM capacity ({Q_max_mvar} MVAR).")
    print("-" * 30)

    # Step 4: Calculate System Losses
    # Post-compensation voltage at Bus B
    V_B_final = cmath.rect(V_B_target_pu, -delta_B)

    # Assuming ideal TR1, V_A_final = V_B_final
    V_A_final = V_B_final
    V_grid = cmath.rect(1.0, 0)
    
    # Current from grid in positive sequence network
    I_grid_pu = (V_grid - V_A_final) / Z1
    
    # Fundamental power loss = Real power supplied by the grid
    P_loss_fundamental_pu = (V_grid * I_grid_pu.conjugate()).real
    P_loss_fundamental_mw = P_loss_fundamental_pu * S_base

    # Total loss including harmonics
    P_loss_total_mw = P_loss_fundamental_mw * harmonic_loss_factor

    print("Step 4: System Power Loss Calculation")
    print(f"Current from grid I_grid = {abs(I_grid_pu):.4f} p.u.")
    print(f"Fundamental power loss = {P_loss_fundamental_mw:.4f} MW")
    print(f"Total power loss (including 4% harmonic losses) = {P_loss_fundamental_mw:.2f} MW * 1.04 = {P_loss_total_mw:.4f} MW")
    print("-" * 30)

    print("\nFinal Results:")
    print(f"The minimum reactive power injection required, Q_opt, is {Q_opt_mvar:.2f} MVAR.")
    print(f"The system's real power losses are {P_loss_total_mw:.2f} MW.")
    
solve_hvac_system()
print(f'<<<52.14>>>')
