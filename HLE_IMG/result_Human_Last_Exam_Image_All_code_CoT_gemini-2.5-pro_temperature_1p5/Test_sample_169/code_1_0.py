import numpy as np
from scipy.optimize import fsolve

def solve_wind_power_optimization():
    """
    Formulates and solves the nonlinear optimization problem for reactive power compensation.
    """

    # 1. System Modeling & Parameter Conversion
    # Base values
    S_base = 100e6  # VA
    V_base = 400e3  # V (at Bus-W)
    Z_base = V_base**2 / S_base  # Ohms

    # Given parameters (assuming they are in p.u., as ohmic values lead to non-physical results)
    R_line_pu = 0.01
    X_line_pu = 0.05
    R_fault_pu = 0.1
    
    # Harmonic losses increase resistance by 6%
    loss_factor = 1.06
    R_line_eff_pu = R_line_pu * loss_factor
    # Ambiguity exists, but assume fault resistance is also part of system losses
    R_fault_eff_pu = R_fault_pu * loss_factor

    # Line impedance from Bus-W to Fault (Z_WF)
    # Assumption: Fault is at midpoint (x=0.5), so Z_GF = Z_WF
    Z_WF_pu = complex(R_line_eff_pu, X_line_pu)
    Z_GF_pu = Z_WF_pu

    # Target voltage and Grid voltage (assumed as reference)
    V_W_target_pu = 0.575
    V_G_pu = 1.0

    # Constraints
    PF_min = 0.95
    tan_phi_max = np.tan(np.arccos(PF_min))
    Q_max_MVAR = 10
    Q_max_pu = Q_max_MVAR / (S_base / 1e6)

    # 2. Network Analysis (Y-bus method)
    y_WF_pu = 1 / Z_WF_pu
    y_GF_pu = 1 / Z_GF_pu
    y_F0_pu = 1 / R_fault_eff_pu

    # Y-bus elements for the W-F-G network
    # From I_F = Y_FW*V_W + Y_FF*V_F + Y_FG*V_G = 0 (no injection at F)
    # V_F = -(Y_FW*V_W + Y_FG*V_G) / Y_FF
    # Y_FW = -y_WF, Y_FG = -y_GF, Y_FF = y_WF + y_GF + y_F0
    Y_FF = y_WF_pu + y_GF_pu + y_F0_pu
    
    # Power injection at Bus W: S_inj_W = V_W * I_W*
    # I_W = Y_WW*V_W + Y_WF*V_F = y_WF*V_W - y_WF*V_F
    # After substituting V_F: S_inj_W = Y_A*|V_W|^2 - Y_B*V_W*V_G*
    # where Y_A = y_WF* (1 + y_WF/Y_FF) and Y_B = y_WF*y_GF/Y_FF
    # Note: These Y_A, Y_B are different from manual calculation steps,
    # let's use the explicit power equations derived from power flow for clarity.
    
    V_W_phasor = V_W_target_pu # Magnitude, angle is unknown 'delta'
    V_G_phasor = complex(V_G_pu, 0) # Reference angle

    def get_power(delta):
        V_W = V_W_phasor * (np.cos(delta) + 1j * np.sin(delta))
        V_F = (y_WF_pu * V_W + y_GF_pu * V_G_phasor) / Y_FF
        I_W = y_WF_pu * (V_W - V_F)
        S_W_inj = V_W * I_W.conjugate()
        return S_W_inj.real, S_W_inj.imag

    # 3. Formulate and Solve Optimization
    # Constraint: Q_eff / P_W <= tan_phi_max  => tan_phi_max * P_W - Q_eff >= 0
    def pf_constraint(delta):
        P_W, Q_eff = get_power(delta)
        return tan_phi_max * P_W - Q_eff

    # The optimal solution under a "minimize Q" objective often lies on the boundary
    # of the feasible region. Here, the binding constraint is the power factor limit.
    # We find the delta values where the PF constraint is met with equality.
    
    # Solve pf_constraint(delta) = 0 to find the boundary angles
    # Initial guesses for delta (usually a small negative angle and a larger positive angle)
    try:
        delta_sol1 = fsolve(pf_constraint, -0.5)[0]
        delta_sol2 = fsolve(pf_constraint, 1.0)[0]
    except Exception:
        print("Could not find a solution for the power factor constraint boundary.")
        return

    # 4. Evaluate solutions and select the optimal one
    P1, Q1 = get_power(delta_sol1)
    P2, Q2 = get_power(delta_sol2)

    # We need to inject reactive power, so Q_eff must be positive.
    # We also want to minimize the injection.
    potential_solutions = []
    if Q1 >= 0:
        potential_solutions.append(Q1)
    if Q2 >= 0:
        potential_solutions.append(Q2)

    if not potential_solutions:
        print("No solution found that provides positive reactive power injection while meeting all constraints.")
        print("The system may be unable to reach the target voltage and power factor simultaneously.")
        return

    # Optimal Q is the minimum valid solution
    Q_opt_pu = min(potential_solutions)

    # Check against compensator capacity
    if Q_opt_pu > Q_max_pu:
        print(f"Optimal reactive power required is {Q_opt_pu:.4f} p.u. ({Q_opt_pu*100:.2f} MVAR).")
        print(f"This exceeds the maximum capacity of the compensator ({Q_max_pu:.4f} p.u. or {Q_max_MVAR} MVAR).")
        print("Therefore, there is no feasible solution with the given compensator.")
        # As per the task, we should still output the calculated optimal value.
    else:
        print(f"Optimal reactive power injection found: {Q_opt_pu:.4f} p.u. ({Q_opt_pu*100:.2f} MVAR).")
        print("This value is within the compensator's capacity.")

    # Find the P and V that correspond to the optimal Q
    final_delta = delta_sol1 if Q1 == Q_opt_pu else delta_sol2
    P_final_pu, Q_final_pu = get_power(final_delta)
    
    # Print the final results in equation form as requested
    V_final_kV = V_W_target_pu * V_base / 1000
    P_final_MW = P_final_pu * S_base / 1e6
    Q_final_MVAR = Q_final_pu * S_base / 1e6

    print("\n--- Final Optimal State ---")
    print(f"Voltage at Bus-W (V_W) = {V_W_target_pu:.3f} p.u. = {V_final_kV:.1f} kV")
    print(f"Active Power from Wind Farm (P_W) = {P_final_pu:.4f} p.u. = {P_final_MW:.2f} MW")
    print(f"Optimal Reactive Power Injection (Q_opt) = {Q_final_pu:.4f} p.u. = {Q_final_MVAR:.2f} MVAR")

    final_S_MVA = np.sqrt(P_final_MW**2 + Q_final_MVAR**2)
    final_PF = P_final_MW / final_S_MVA
    
    print(f"Resulting Power Factor = {final_PF:.4f} lagging")
    
    print("\nFinal equation based on the optimization result:")
    print(f"To restore V_W to {V_W_target_pu} kV while maintaining PF >= {PF_min},")
    print("the required reactive power injection is:")
    print(f"Q_opt = {Q_final_MVAR:.4f} MVAR")
    
    # Final answer format
    print(f"\n<<< {Q_final_MVAR:.4f} >>>")


solve_wind_power_optimization()