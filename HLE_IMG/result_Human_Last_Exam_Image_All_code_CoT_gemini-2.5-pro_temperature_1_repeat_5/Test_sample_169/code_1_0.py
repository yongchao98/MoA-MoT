import cmath
import math

def solve_wind_compensation():
    """
    Calculates the optimal reactive power injection to mitigate a fault in a wind energy system.
    """
    # 1. Define system parameters in per-unit (p.u.)
    S_base_MVA = 100.0  # MVA
    V_W_target_pu = 0.575  # p.u.
    
    R_WF_pu = 0.01  # p.u.
    X_WF_pu = 0.05  # p.u.
    Z_WF_pu = complex(R_WF_pu, X_WF_pu)
    
    R_F_pu = 0.1  # p.u.
    
    Q_max_MVAR = 10.0
    Q_max_pu = Q_max_MVAR / S_base_MVA
    
    PF_min = 0.95  # lagging
    harmonic_loss_increase = 0.06 # 6%

    print("--- Step 1: System Parameters ---")
    print(f"Base Power (S_base): {S_base_MVA} MVA")
    print(f"Target Voltage at Bus-W (V_W): {V_W_target_pu} p.u.")
    print(f"Transmission Impedance (Z_WF): {R_WF_pu} + j{X_WF_pu} p.u.")
    print(f"Fault Resistance (R_F): {R_F_pu} p.u.")
    print(f"Harmonic Loss Factor: 1 + {harmonic_loss_increase} = {1 + harmonic_loss_increase}")
    print(f"Minimum Power Factor (PF_min): {PF_min} lagging")
    print(f"Compensator Max Capacity (Q_max): {Q_max_MVAR} MVAR ({Q_max_pu} p.u.)")
    print("-" * 35 + "\n")

    # 2. Calculate effective impedance including harmonic losses
    R_total_pu = R_WF_pu + R_F_pu
    R_eff_pu = R_total_pu * (1 + harmonic_loss_increase)
    Z_eff_pu = complex(R_eff_pu, X_WF_pu)

    print("--- Step 2: Calculate Effective Impedance ---")
    print(f"Total fundamental resistance R_total = R_WF + R_F = {R_WF_pu} + {R_F_pu} = {R_total_pu:.4f} p.u.")
    print(f"Effective resistance with harmonic losses R_eff = R_total * 1.06 = {R_total_pu:.4f} * {1 + harmonic_loss_increase} = {R_eff_pu:.4f} p.u.")
    print(f"Effective fault impedance Z_eff = {Z_eff_pu.real:.4f} + j{Z_eff_pu.imag:.4f} p.u.")
    print("-" * 35 + "\n")

    # 3. Calculate required power S_W = P_W + jQ_W at Bus-W
    # S_W = |V_W|^2 / conj(Z_eff)
    V_W_mag_sq = V_W_target_pu ** 2
    S_W_pu = V_W_mag_sq / Z_eff_pu.conjugate()
    P_W_pu = S_W_pu.real
    Q_W_pu = S_W_pu.imag

    print("--- Step 3: Calculate Required Power at Bus-W ---")
    print("This is the power drawn by the faulted line to maintain V_W = 0.575 p.u.")
    print(f"S_W = |V_W|^2 / Z_eff* = {V_W_target_pu}^2 / ({Z_eff_pu.real:.4f} - j{Z_eff_pu.imag:.4f})")
    print(f"S_W = {P_W_pu:.4f} + j{Q_W_pu:.4f} p.u.")
    print(f"Required Real Power (P_W): {P_W_pu:.4f} p.u.")
    print(f"Required Reactive Power (Q_W): {Q_W_pu:.4f} p.u.")
    print("-" * 35 + "\n")
    
    # 4. Formulate and solve the optimization
    # To minimize Q_comp, we must maximize Q_gen subject to the PF constraint
    # PF_gen >= 0.95  =>  |Q_gen / P_gen| <= tan(acos(0.95))
    
    tan_phi_max = math.tan(math.acos(PF_min))
    P_gen_pu = P_W_pu # Generator must supply the real power
    Q_gen_max_pu = P_gen_pu * tan_phi_max
    
    # Q_opt is the minimum Q_comp needed
    # Q_W = Q_gen + Q_comp  => Q_comp = Q_W - Q_gen
    # To minimize Q_comp, use Q_gen_max
    Q_opt_pu = Q_W_pu - Q_gen_max_pu

    print("--- Step 4: Determine Optimal Reactive Power Injection ---")
    print(f"The generator must supply all the real power: P_gen = P_W = {P_gen_pu:.4f} p.u.")
    print(f"From the PF constraint (PF >= {PF_min}), the max generator reactive power is:")
    print(f"Q_gen_max = P_gen * tan(acos({PF_min})) = {P_gen_pu:.4f} * {tan_phi_max:.4f} = {Q_gen_max_pu:.4f} p.u.")
    print("\nThe optimal (minimum) compensator injection Q_opt is the total required Q_W minus what the generator can supply:")
    print(f"Q_opt = Q_W - Q_gen_max")
    print(f"Q_opt = {Q_W_pu:.4f} - {Q_gen_max_pu:.4f} = {Q_opt_pu:.4f} p.u.")
    
    Q_opt_MVAR = Q_opt_pu * S_base_MVA
    print(f"\nIn physical units, Q_opt = {Q_opt_pu:.4f} p.u. * {S_base_MVA} MVA = {Q_opt_MVAR:.2f} MVAR.")
    print("-" * 35 + "\n")

    # 5. Final Feasibility Check
    print("--- Step 5: Feasibility Check ---")
    print(f"The required optimal injection is Q_opt = {Q_opt_MVAR:.2f} MVAR.")
    print(f"The maximum capacity of the compensating device is Q_max = {Q_max_MVAR:.2f} MVAR.")
    
    if Q_opt_MVAR <= Q_max_MVAR:
        print("\nConclusion: The solution is feasible. The device can supply the required reactive power.")
    else:
        print("\nConclusion: The solution is not feasible. The required reactive power ({:.2f} MVAR) exceeds the device's maximum capacity ({:.2f} MVAR).".format(Q_opt_MVAR, Q_max_MVAR))
    print("-" * 35 + "\n")
    
    print(f"Final determined optimal reactive power injection: {Q_opt_MVAR:.2f} MVAR.")

solve_wind_compensation()