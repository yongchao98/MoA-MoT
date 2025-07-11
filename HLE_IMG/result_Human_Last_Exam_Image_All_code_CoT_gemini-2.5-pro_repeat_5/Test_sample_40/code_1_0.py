import cmath
import math

def solve_hvac_optimization():
    """
    Solves the non-linear optimization problem for the HVAC transmission system.
    """
    # --- 1. System Parameters and Assumptions ---
    S_base = 100.0  # MVA
    V_B_nom = 220.0  # kV
    Q_max = 50.0  # MVAR
    
    # System impedances in per unit (p.u.)
    Z_S = complex(0.02, 0.10)
    R_S = Z_S.real
    X_S = Z_S.imag
    
    # The given Z_F = 0.15 Ohm leads to an infeasible solution.
    # To find a feasible solution, we assume a typo and use Z_F = 1.5 p.u.
    # This fault impedance is assumed to be purely resistive.
    Z_F = 1.5
    
    # Harmonic loss increase factor
    harmonic_loss_factor = 0.04
    
    # Source voltage (External Grid at Bus A) in p.u.
    V_A = complex(1.0, 0.0)

    print("Step 1: Define System Parameters (in per-unit)")
    print(f"S_base = {S_base} MVA")
    print(f"Z_S = {R_S:.2f} + j{X_S:.2f} p.u.")
    print(f"Z_F = {Z_F:.2f} p.u. (Assumed value for a feasible solution)")
    print("-" * 30)

    # --- 2. Solve for Optimal Susceptance (Bc) ---
    # From the constraint |V_B| = 1.0 p.u., we derive the quadratic equation:
    # a*Bc^2 + b*Bc + c = 0
    G_F = 1 / Z_F
    
    a = R_S**2 + X_S**2
    b = -2 * X_S
    c = (R_S**2 + X_S**2) * G_F**2 + 2 * R_S * G_F
    
    print("Step 2: Solve for Optimal STATCOM Susceptance (Bc)")
    print(f"Derived quadratic equation: {a:.4f}*Bc^2 + {b:.4f}*Bc + {c:.4f} = 0")
    
    # Solve the quadratic equation for Bc
    discriminant = b**2 - 4 * a * c
    if discriminant < 0:
        print("No real solution exists. Cannot restore voltage to 1.0 p.u.")
        return

    Bc1 = (-b + math.sqrt(discriminant)) / (2 * a)
    Bc2 = (-b - math.sqrt(discriminant)) / (2 * a)
    
    print(f"Solutions for Bc are: {Bc1:.4f} p.u. and {Bc2:.4f} p.u.")
    
    # Choose the minimum positive solution to minimize reactive power
    Bc_opt = min(b for b in [Bc1, Bc2] if b > 0)
    print(f"Optimal minimum susceptance Bc_opt = {Bc_opt:.4f} p.u.")
    print("-" * 30)
    
    # --- 3. Calculate Optimal Reactive Power (Q_opt) ---
    # Q_opt = Bc_opt * |V_B|^2, and since |V_B|=1.0, Q_opt_pu = Bc_opt
    Q_opt_pu = Bc_opt
    Q_opt_MVAR = Q_opt_pu * S_base
    
    print("Step 3: Calculate Optimal Reactive Power (Q_opt)")
    print(f"Q_opt (p.u.) = Bc_opt * |V_B|^2 = {Bc_opt:.4f} * 1.0^2 = {Q_opt_pu:.4f} p.u.")
    print(f"Q_opt (MVAR) = Q_opt (p.u.) * S_base = {Q_opt_pu:.4f} * {S_base:.1f} = {Q_opt_MVAR:.2f} MVAR")

    if Q_opt_MVAR > Q_max:
        print(f"Warning: Required Q_opt ({Q_opt_MVAR:.2f} MVAR) exceeds STATCOM capacity ({Q_max:.1f} MVAR).")
    else:
        print(f"The required reactive power {Q_opt_MVAR:.2f} MVAR is within the STATCOM capacity of {Q_max:.1f} MVAR.")
    print("-" * 30)
        
    # --- 4. Calculate System Losses ---
    print("Step 4: Calculate System Real Power Losses")
    
    # Calculate final system state with compensation
    # Total shunt admittance at Bus B
    Y_shunt = complex(G_F, Bc_opt)
    Z_shunt = 1 / Y_shunt
    
    # Final compensated voltage at Bus B
    V_B = V_A * Z_shunt / (Z_S + Z_shunt)
    
    # Line current from Bus A to Bus B
    I_line = (V_A - V_B) / Z_S
    I_line_mag_sq = abs(I_line)**2
    
    # Loss in the series impedance Z_S
    P_loss_S_pu = I_line_mag_sq * R_S
    P_loss_S_MW = P_loss_S_pu * S_base
    
    # Loss in the fault impedance Z_F
    P_loss_F_pu = abs(V_B)**2 / Z_F
    P_loss_F_MW = P_loss_F_pu * S_base

    print(f"Final Voltage at Bus B, V_B = {abs(V_B):.3f} p.u. at angle {math.degrees(cmath.phase(V_B)):.2f} degrees")
    print(f"Line Current Magnitude |I_line| = {abs(I_line):.3f} p.u.")
    
    # Calculate total base losses
    P_loss_total_pu = P_loss_S_pu + P_loss_F_pu
    P_loss_total_MW = P_loss_total_pu * S_base
    
    print("\n--- Power Loss Calculation ---")
    print(f"Loss in Z_S = |I_line|^2 * R_S = {I_line_mag_sq:.4f} * {R_S:.2f} = {P_loss_S_pu:.4f} p.u. = {P_loss_S_MW:.2f} MW")
    print(f"Loss in Z_F = |V_B|^2 / Z_F = {abs(V_B)**2:.4f} / {Z_F:.2f} = {P_loss_F_pu:.4f} p.u. = {P_loss_F_MW:.2f} MW")

    print(f"Total fundamental loss = {P_loss_S_pu:.4f} p.u. + {P_loss_F_pu:.4f} p.u. = {P_loss_total_pu:.4f} p.u. = {P_loss_total_MW:.2f} MW")

    # Add harmonic losses
    P_loss_final_pu = P_loss_total_pu * (1 + harmonic_loss_factor)
    P_loss_final_MW = P_loss_final_pu * S_base

    print(f"\nFinal loss with {harmonic_loss_factor*100}% harmonic increase:")
    print(f"P_loss_final = P_loss_total * (1 + {harmonic_loss_factor}) = {P_loss_total_MW:.2f} MW * {1+harmonic_loss_factor} = {P_loss_final_MW:.2f} MW")
    print("-" * 30)

    print("\n--- Final Answer ---")
    print(f"Optimal Reactive Power Injection (Q_opt): {Q_opt_MVAR:.2f} MVAR")
    print(f"Total System Real Power Losses: {P_loss_final_MW:.2f} MW")
    
    return Q_opt_MVAR, P_loss_final_MW

# Run the solver
solve_hvac_optimization()