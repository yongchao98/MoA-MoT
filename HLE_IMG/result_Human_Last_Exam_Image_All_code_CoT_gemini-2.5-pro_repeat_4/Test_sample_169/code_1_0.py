import math

def solve_reactive_power_optimization():
    """
    Solves the nonlinear optimization problem to find the optimal reactive power injection.
    """
    # Step 1: Define system parameters and establish the per-unit system
    S_base = 100.0  # MVA
    V_base = 0.575  # kV (chosen to match the target voltage at Bus-W)
    V_W_target_kV = 0.575 # kV
    Z_base = (V_base * 1000)**2 / (S_base * 10**6)  # in Ohms

    Z_WF_ohm = complex(0.01, 0.05)  # Ohms
    R_F_ohm = 0.1  # Ohms
    Q_max_MVAR = 10.0  # MVAR
    harmonic_loss_increase = 0.06  # 6%
    min_power_factor = 0.95  # lagging

    # Step 2: Convert to per-unit
    V_W_target_pu = V_W_target_kV / V_base
    Z_WF_pu = Z_WF_ohm / Z_base
    R_F_pu = R_F_ohm / Z_base
    Q_max_pu = Q_max_MVAR / S_base

    # Step 3: Simplified fault model
    Z_total_pu = Z_WF_pu + R_F_pu

    # Step 4: Calculate required fundamental power S_W = P_fund + jQ_net
    # S_W = |V_W|^2 / Z_total_pu* (where * denotes complex conjugate)
    S_W_fund_pu = (V_W_target_pu**2) / Z_total_pu.conjugate()
    P_fund_pu = S_W_fund_pu.real
    Q_net_pu = S_W_fund_pu.imag

    # Step 5: Account for harmonic losses
    P_gen_pu = P_fund_pu * (1 + harmonic_loss_increase)

    # Step 6: Apply power factor constraint to find max Q_gen
    # For a lagging PF, Q_gen >= 0. The constraint is PF = P_gen / |S_gen| >= 0.95
    # This is equivalent to |Q_gen| <= P_gen * tan(acos(PF))
    phi_max = math.acos(min_power_factor)
    tan_phi_max = math.tan(phi_max)
    Q_gen_max_pu = P_gen_pu * tan_phi_max

    # Step 7: Determine optimal Q_comp
    # To minimize Q_comp, we maximize Q_gen.
    # Q_comp = Q_net - Q_gen
    Q_opt_pu = Q_net_pu - Q_gen_max_pu

    # Step 8: Convert back to MVAR and present the result
    Q_opt_MVAR = Q_opt_pu * S_base

    print("--- Calculation Steps ---")
    print(f"1. Base Impedance (Z_base): {Z_base:.5f} Ohms")
    print(f"2. Total impedance to fault (Z_total): {Z_total_pu.real:.3f} + j{Z_total_pu.imag:.3f} p.u.")
    print(f"3. Required fundamental power at Bus-W (S_W_fund): {P_fund_pu:.5f} + j{Q_net_pu:.5f} p.u.")
    print(f"4. Total active power required from generator (P_gen), including harmonic losses: {P_fund_pu:.5f} * (1 + {harmonic_loss_increase}) = {P_gen_pu:.5f} p.u.")
    print(f"5. Max reactive power from generator (Q_gen_max) at PF={min_power_factor}: {P_gen_pu:.5f} * tan(acos({min_power_factor})) = {Q_gen_max_pu:.5f} p.u.")
    
    print("\n--- Final Equation ---")
    print("The optimal reactive power injection (Q_opt) is the total reactive power needed minus the maximum reactive power the generator can supply.")
    print("Q_opt (p.u.) = Q_net_required (p.u.) - Q_gen_max (p.u.)")
    print(f"Q_opt (p.u.) = {Q_net_pu:.5f} - {Q_gen_max_pu:.5f} = {Q_opt_pu:.5f} p.u.")
    
    print("\n--- Final Answer ---")
    print(f"Optimal reactive power injection Q_opt = {Q_opt_pu:.5f} p.u. * {S_base} MVA = {Q_opt_MVAR:.3f} MVAR")
    
    # Check if solution is within device limits
    if Q_opt_MVAR > Q_max_MVAR:
        print(f"\nWarning: Optimal value {Q_opt_MVAR:.3f} MVAR exceeds device max capacity of {Q_max_MVAR} MVAR.")
    else:
        print(f"\nThe calculated value {Q_opt_MVAR:.3f} MVAR is within the device's maximum capacity of {Q_max_MVAR} MVAR.")
    
    # Return the final numerical answer as requested
    final_answer = round(Q_opt_MVAR, 3)
    return final_answer

# Execute the function and print the final answer in the required format
optimal_q = solve_reactive_power_optimization()
print(f"\n<<<{optimal_q}>>>")
