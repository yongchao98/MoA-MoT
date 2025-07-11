import math

def solve_wind_energy_optimization():
    """
    Solves the nonlinear optimization problem for the wind energy system.
    """
    # Step 1: Define system parameters
    S_base = 100e6  # VA
    V_base = 0.575e3  # V, chosen based on problem statement
    V_W_target_kV = 0.575 # kV
    R_WF_ohm = 0.01  # Ohms
    X_WF_ohm = 0.05  # Ohms
    R_F_ohm = 0.1  # Ohms
    Q_max_MVAR = 10.0  # MVAR
    harmonic_loss_increase = 0.06  # 6%

    # Step 2: Calculate base values and convert to per-unit
    Z_base = V_base**2 / S_base
    R_WF_pu = R_WF_ohm / Z_base
    X_WF_pu = X_WF_ohm / Z_base
    R_F_pu = R_F_ohm / Z_base
    Q_max_pu = Q_max_MVAR * 1e6 / S_base
    V_W_target_pu = V_W_target_kV * 1e3 / V_base # This will be 1.0

    # Step 3 & 4: Model Thevenin impedance including harmonic losses
    # Thevenin source is assumed to be V_th = 1.0 at angle 0.
    # Thevenin impedance is the line impedance.
    R_th_eq_pu = R_WF_pu * (1 + harmonic_loss_increase)
    X_th_pu = X_WF_pu
    Z_th_sq_mag = R_th_eq_pu**2 + X_th_pu**2

    print("--- Per-Unit System Parameters ---")
    print(f"Base Impedance (Z_base): {Z_base:.6f} Ohms")
    print(f"Equivalent Thevenin Resistance (R_th_eq): {R_th_eq_pu:.4f} p.u.")
    print(f"Thevenin Reactance (X_th): {X_th_pu:.4f} p.u.")
    print(f"Fault Resistance (R_F): {R_F_pu:.4f} p.u.")
    print(f"Target Bus-W Voltage (|V_W|): {V_W_target_pu:.2f} p.u.")
    print(f"Maximum Reactive Power (Q_max): {Q_max_pu:.2f} p.u.")
    print("-" * 35)

    # Step 5 & 6: Solve for the bus voltage angle delta_W
    # The real power balance equation is: A*cos(delta) + B*sin(delta) = C
    # where A = R_th_eq_pu, B = -X_th_pu
    # and C = R_th_eq_pu + (Z_th_sq_mag / R_F_pu)
    A = R_th_eq_pu
    B = -X_th_pu
    C = R_th_eq_pu + (Z_th_sq_mag / R_F_pu)

    # Solve R*cos(delta - alpha) = C, where R=sqrt(A^2+B^2) and alpha=atan2(B,A)
    R = math.sqrt(A**2 + B**2)
    alpha = math.atan2(B, A)

    cos_val = C / R
    if abs(cos_val) > 1:
        print("No real solution exists for the voltage angle.")
        return

    theta = math.acos(cos_val)

    # Two possible solutions for delta_W
    delta_W_rad_1 = theta + alpha
    delta_W_rad_2 = -theta + alpha

    solutions = []
    for delta_W in [delta_W_rad_1, delta_W_rad_2]:
        # Step 7: Calculate Q_comp for the angle
        # The imaginary power balance gives: -Q_comp = [X_th*(cos(d)-1) + R_th*sin(d)] / Z_th_sq_mag
        cos_d = math.cos(delta_W)
        sin_d = math.sin(delta_W)
        
        num_q = X_th_pu * (cos_d - 1) + R_th_eq_pu * sin_d
        Q_comp_pu = -num_q / Z_th_sq_mag
        Q_comp_MVAR = Q_comp_pu * (S_base / 1e6)
        
        # Step 8: Check if solution is feasible
        if 0 <= Q_comp_MVAR <= Q_max_MVAR:
            solutions.append({
                "delta_W_deg": math.degrees(delta_W),
                "Q_comp_MVAR": Q_comp_MVAR,
                "Q_comp_pu": Q_comp_pu,
                "cos_d": cos_d,
                "sin_d": sin_d
            })

    # Step 9: Select the optimal solution
    if not solutions:
        print("No feasible solution found that satisfies the constraints.")
    else:
        # The optimal solution is the one that minimizes Q_comp
        optimal_solution = min(solutions, key=lambda x: x["Q_comp_MVAR"])
        Q_opt = optimal_solution["Q_comp_MVAR"]
        
        print("\n--- Optimal Solution Found ---")
        print(f"Bus-W Voltage Angle (delta_W): {optimal_solution['delta_W_deg']:.2f} degrees")
        print(f"Optimal Reactive Power Injection (Q_opt): {Q_opt:.3f} MVAR")

        # Print the final equation with numbers
        print("\n--- Final Calculation Breakdown ---")
        print("The optimal reactive power Q_opt (in p.u.) is calculated using:")
        print("Q_comp_pu = -[X_th * (cos(delta_W) - 1) + R_th_eq * sin(delta_W)] / (R_th_eq^2 + X_th^2)")
        print("Substituting the values:")
        print(f"Q_comp_pu = -[{X_th_pu:.4f} * ({optimal_solution['cos_d']:.4f} - 1) + {R_th_eq_pu:.4f} * ({optimal_solution['sin_d']:.4f})] / {Z_th_sq_mag:.4f}")
        
        num_val = X_th_pu * (optimal_solution['cos_d'] - 1) + R_th_eq_pu * optimal_solution['sin_d']
        final_q_pu = -num_val / Z_th_sq_mag
        
        print(f"Q_comp_pu = -[{num_val:.4f}] / {Z_th_sq_mag:.4f} = {final_q_pu:.5f} p.u.")
        print(f"Q_opt (MVAR) = {final_q_pu:.5f} p.u. * {S_base/1e6:.0f} MVA = {Q_opt:.3f} MVAR")

        print(f"\nThe optimal reactive power injection required is {Q_opt:.3f} MVAR.")
        print(f"<<<{Q_opt:.3f}>>>")

solve_wind_energy_optimization()