import numpy as np

def solve_reactive_power_optimization():
    """
    Formulates and solves the nonlinear optimization problem to find the optimal
    reactive power injection from a fault compensating device in a wind energy system.
    """
    # 1. Define System Parameters and Assumptions
    # Base and target values
    S_base_MVA = 100.0
    V_W_target_pu = 0.575  # Target voltage at Bus-W in per-unit
    V_G = 1.0 + 0j         # Grid voltage at Bus-G (infinite bus) in p.u.

    # Impedances and fault resistance in per-unit (p.u.)
    Z_WF = 0.01 + 0.05j
    Z_GF = 0.01 + 0.05j  # Assumed midpoint fault
    R_F = 0.1

    # Constraints
    Q_max_MVAR = 10.0
    Q_max_pu = Q_max_MVAR / S_base_MVA
    power_factor_min = 0.95
    # The constraint |Q_gen| / P_gen <= tan(acos(PF))
    k_pf = np.tan(np.arccos(power_factor_min))
    
    # Harmonic loss factor
    harmonic_loss_factor = 1.06 # 6% increase

    # 2. Optimization by searching over the voltage angle delta_W
    min_q_comp = float('inf')
    optimal_params = {}

    # Search over a plausible range of angles for power export
    delta_w_range = np.linspace(-np.pi / 2, 0, 10000)

    for delta_w in delta_w_range:
        V_W = V_W_target_pu * np.exp(1j * delta_w)

        # Calculate voltage at the fault point F using KCL
        V_F = (V_W + V_G) * R_F / (2 * R_F + Z_WF)

        # Calculate power flowing out of Bus-W (S_out = P_out + j*Q_out)
        I_WF = (V_W - V_F) / Z_WF
        S_out = V_W * np.conj(I_WF)
        P_out = S_out.real
        Q_out = S_out.imag

        if P_out < 0:
            continue  # We are interested in power export scenarios

        # Determine the required Q_comp to satisfy the generator's PF limit.
        # To minimize Q_comp, the generator should provide maximum possible reactive power,
        # so we set the generator's reactive power Q_W to its limit: Q_W = k_pf * P_out.
        # Since Q_out = Q_W + Q_comp, the required Q_comp is Q_out - Q_W.
        q_comp_required = Q_out - k_pf * P_out
        
        # The compensator cannot absorb reactive power, so Q_comp must be >= 0.
        q_comp_chosen = max(0, q_comp_required)

        # Check if this solution is feasible within the compensator's capacity
        if q_comp_chosen <= Q_max_pu:
            if q_comp_chosen < min_q_comp:
                min_q_comp = q_comp_chosen
                optimal_params = {
                    'delta_w_deg': np.rad2deg(delta_w),
                    'Q_comp_pu': q_comp_chosen,
                    'P_W_pu': P_out,
                    'Q_out_pu': Q_out,
                    'V_F': V_F,
                }

    # 3. Print the results
    if not optimal_params:
        print("No feasible solution found that satisfies all constraints.")
        return

    print("--- Optimal Solution Found ---")
    Q_opt_pu = optimal_params['Q_comp_pu']
    Q_opt_MVAR = Q_opt_pu * S_base_MVA
    P_W_pu = optimal_params['P_W_pu']
    Q_out_pu = optimal_params['Q_out_pu']
    V_F_opt = optimal_params['V_F']
    
    # Calculate generator's reactive power at the optimal point
    Q_W_pu = Q_out_pu - Q_opt_pu
    
    # Calculate final power factor of the generator
    final_pf = P_W_pu / np.sqrt(P_W_pu**2 + Q_W_pu**2)

    # Calculate system losses
    I_GF = (V_F_opt - V_G) / Z_GF
    P_G_out_pu = - (V_G * np.conj(I_GF)).real # Power delivered to grid
    
    fundamental_loss_pu = P_W_pu - P_G_out_pu
    total_loss_pu = fundamental_loss_pu * harmonic_loss_factor
    total_loss_MW = total_loss_pu * S_base_MVA
    
    print(f"\nOptimal reactive power injection Q_opt = {Q_opt_MVAR:.2f} MVAR ({Q_opt_pu:.4f} p.u.)")
    print("\nThis result is determined by solving the equation:")
    print("Q_opt = max(0, Q_out - k_pf * P_W)")
    print("\nAt the optimal point:")
    print(f"  Total real power from wind farm (P_W) = {P_W_pu:.4f} p.u.")
    print(f"  Total reactive power from Bus-W (Q_out) = {Q_out_pu:.4f} p.u.")
    print(f"  Power factor constant (k_pf = tan(acos({power_factor_min}))) = {k_pf:.4f}")
    
    print("\nSubstituting these values into the equation:")
    print(f"  Q_opt = max(0, {Q_out_pu:.4f} - {k_pf:.4f} * {P_W_pu:.4f})")
    print(f"  Q_opt = max(0, {Q_out_pu:.4f} - {k_pf * P_W_pu:.4f})")
    print(f"  Q_opt = {Q_opt_pu:.4f} p.u.")

    print("\n--- System State at Optimal Point ---")
    print(f"Wind Farm Real Power (P_W): {P_W_pu * S_base_MVA:.2f} MW")
    print(f"Wind Farm Reactive Power (Q_W): {Q_W_pu * S_base_MVA:.2f} MVAR")
    print(f"Wind Farm Power Factor: {final_pf:.4f}")
    print(f"Total System Losses (incl. harmonics): {total_loss_MW:.2f} MW")
    
if __name__ == '__main__':
    solve_reactive_power_optimization()
