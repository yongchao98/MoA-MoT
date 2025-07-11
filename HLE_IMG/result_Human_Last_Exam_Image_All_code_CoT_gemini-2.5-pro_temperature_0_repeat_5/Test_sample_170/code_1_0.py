import cmath
import math

def solve_hvac_optimization():
    """
    Solves the non-linear optimization problem for the HVAC transmission system.
    """
    # --- Given Data and Base Values ---
    S_base_MVA = 100.0  # Base power in MVA
    V_B_nom_kV = 220.0  # Nominal voltage at Bus B in kV
    
    # System impedances (Assumed to be in per-unit, as Ohm values are not feasible)
    Z_S = 0.02 + 0.10j  # Transmission line impedance in pu
    Z_F = 0.15 + 0.0j   # Fault impedance in pu (purely resistive)
    
    # System conditions
    V_A_pu = 1.0             # Voltage at Bus A (source) in pu
    V_B_fault_pu = 0.85      # Voltage at Bus B during fault in pu
    V_B_target_pu = 1.0      # Target voltage at Bus B after compensation in pu
    loss_increase_factor = 1.04 # 4% increase in losses due to harmonics
    
    # STATCOM constraints
    Q_max_MVAR = 50.0
    
    print("--- Step 1: Estimate Fault Location 'x' ---")
    # We solve for x in |V_B| = |V_A * Z_F / (x*Z_S + Z_F)|
    # This simplifies to a quadratic equation for x: a*x^2 + b*x + c = 0
    a = (abs(Z_S)**2) * (V_B_fault_pu**2)
    b = 2 * (Z_S.real * Z_F.real + Z_S.imag * Z_F.imag) * (V_B_fault_pu**2)
    c = (abs(Z_F)**2) * (V_B_fault_pu**2) - (abs(V_A_pu * Z_F)**2)
    
    # Solve the quadratic equation for x
    discriminant = b**2 - 4*a*c
    if discriminant < 0:
        print("Error: Cannot find a real fault location.")
        return
        
    x1 = (-b + math.sqrt(discriminant)) / (2*a)
    x2 = (-b - math.sqrt(discriminant)) / (2*a)
    
    # Choose the physically meaningful solution (0 < x < 1)
    x = x1 if 0 < x1 < 1 else x2
    print(f"Estimated fault location x = {x:.4f} (fraction of line length from Bus A)")

    print("\n--- Step 2: Calculate Thevenin Equivalent at Bus B ---")
    # Z_th = (1-x)*Z_S + (x*Z_S * Z_F) / (x*Z_S + Z_F)
    Z_th_parallel = (x * Z_S * Z_F) / (x * Z_S + Z_F)
    Z_th = (1 - x) * Z_S + Z_th_parallel
    R_th = Z_th.real
    X_th = Z_th.imag
    V_th_pu = V_B_fault_pu # Thevenin voltage is the open-circuit voltage at Bus B
    
    print(f"Thevenin Impedance Z_th = {R_th:.4f} + j{X_th:.4f} pu")
    print(f"Thevenin Voltage V_th = {V_th_pu:.4f} pu")

    print("\n--- Step 3: Calculate Optimal Reactive Power Q_opt ---")
    # Solve the quadratic equation for Q_opt: a*Q^2 + b*Q + c = 0
    # from |V_B|^2 = |V_th|^2 + 2*Q*X_th - |Z_th|^2 * Q^2 / |V_B|^2
    a_q = abs(Z_th)**2 / (V_B_target_pu**2)
    b_q = -2 * X_th
    c_q = V_th_pu**2 - V_B_target_pu**2
    
    # Solve for Q_opt
    discriminant_q = b_q**2 - 4*a_q*c_q
    if discriminant_q < 0:
        print("Error: Voltage cannot be restored to the target value.")
        return
        
    Q_sol1 = (-b_q + math.sqrt(discriminant_q)) / (2*a_q)
    Q_sol2 = (-b_q - math.sqrt(discriminant_q)) / (2*a_q)
    
    # The minimum required injection corresponds to the smaller, stable solution
    Q_opt_pu = min(Q_sol1, Q_sol2)
    Q_opt_MVAR = Q_opt_pu * S_base_MVA
    
    print(f"Equation for Q_opt: {a_q:.5f}*Q^2 {b_q:.5f}*Q + {c_q:.5f} = 0")
    print(f"Solutions for Q_opt are {Q_sol1:.4f} pu and {Q_sol2:.4f} pu")
    print(f"The minimum required reactive power is Q_opt = {Q_opt_pu:.4f} pu")
    print(f"This corresponds to Q_opt = {Q_opt_MVAR:.2f} MVAR")
    if Q_opt_MVAR > Q_max_MVAR:
        print(f"Note: This required power ({Q_opt_MVAR:.2f} MVAR) exceeds the STATCOM's maximum capacity of {Q_max_MVAR} MVAR.")

    print("\n--- Step 4: Calculate System Real Power Losses ---")
    # Base loss = |I|^2 * R_th, where I = Q_opt / V_B_target
    current_pu = Q_opt_pu / V_B_target_pu
    P_loss_base_pu = (current_pu**2) * R_th
    
    # Final loss with 4% increase for harmonics
    P_loss_final_pu = P_loss_base_pu * loss_increase_factor
    P_loss_final_MW = P_loss_final_pu * S_base_MVA
    
    print(f"Base system losses = {P_loss_base_pu:.4f} pu")
    print(f"Final system losses (with 4% harmonic increase) = {P_loss_final_pu:.4f} pu")
    print(f"This corresponds to final losses = {P_loss_final_MW:.2f} MW")

    print("\n--- Final Answer ---")
    print(f"The optimal reactive power injection required is Q_opt = {Q_opt_MVAR:.2f} MVAR.")
    print(f"The system's real power losses under these conditions are {P_loss_final_MW:.2f} MW.")

# Execute the function
solve_hvac_optimization()