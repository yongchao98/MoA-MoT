import numpy as np
from scipy.optimize import minimize

def solve_fault_compensation():
    """
    Formulates and solves the nonlinear optimization problem for wind farm fault compensation.
    """
    # --- Step 1: Define System Parameters in Per-Unit (p.u.) ---
    S_base = 100.0      # Base Power in MVA
    V_W_mag = 0.575     # Target voltage magnitude at Bus-W [p.u.]
    V_G_mag = 1.0       # Grid voltage magnitude at Bus-G [p.u.]
    Z_WF = 0.01 + 0.05j # Impedance from Bus-W to Fault [p.u.]
    Z_GF = 0.01 + 0.05j # Assumed: Impedance from Fault to Bus-G [p.u.]
    R_F = 0.1           # Fault resistance [p.u.]
    
    # Convert impedances to admittances
    Y_WF = 1.0 / Z_WF
    Y_GF = 1.0 / Z_GF
    Y_F = 1.0 / R_F
    # Total admittance for nodal analysis at the fault point
    Y_sum_node_F = Y_WF + Y_GF + Y_F

    # Define optimization constraints
    Q_max_pu = 10.0 / S_base  # Max reactive power injection [p.u.]
    PF_min = 0.95             # Minimum power factor (lagging)
    tan_phi_max = np.tan(np.arccos(PF_min)) # |Q/P| <= tan(phi)

    # --- Step 2: Define Functions for Power Calculation ---
    def get_system_state(delta_W):
        """Calculates key system values based on the Bus-W voltage angle."""
        V_W = V_W_mag * (np.cos(delta_W) + 1j * np.sin(delta_W))
        V_G = V_G_mag + 0j

        # Using nodal analysis to find the voltage at the fault point F
        V_F = (V_W * Y_WF + V_G * Y_GF) / Y_sum_node_F
        
        # Calculate complex power injected at Bus-W
        I_W_inj = (V_W - V_F) * Y_WF
        S_W_inj = V_W * np.conj(I_W_inj)
        
        return S_W_inj, V_W, V_F, V_G

    # --- Step 3: Define and Solve the Optimization Problem ---
    # Objective: Minimize reactive power Q_comp = S_W_inj.imag
    def objective_func(delta_W):
        S_W, _, _, _ = get_system_state(delta_W)
        return S_W.imag

    # Constraints (must be >= 0 for SLSQP method)
    constraints = [
        # 1. Q_comp >= 0
        {'type': 'ineq', 'fun': lambda d: objective_func(d)},
        # 2. Q_comp <= Q_max_pu
        {'type': 'ineq', 'fun': lambda d: Q_max_pu - objective_func(d)},
        # 3. P_W >= 0
        {'type': 'ineq', 'fun': lambda d: get_system_state(d)[0].real},
        # 4. PF >= 0.95 lagging (Q/P <= tan_phi_max)
        {'type': 'ineq', 'fun': lambda d: tan_phi_max * get_system_state(d)[0].real - objective_func(d)}
    ]
    
    # Initial guess and bounds for the decision variable delta_W
    initial_delta_W = -0.5
    bounds = [(-np.pi, np.pi)]

    solution = minimize(objective_func, initial_delta_W, method='SLSQP', bounds=bounds, constraints=constraints)

    # --- Step 4: Output and Analyze Results ---
    if solution.success:
        Q_opt_pu = solution.fun
        Q_opt_mvar = Q_opt_pu * S_base
        
        optimal_delta = solution.x[0]
        S_W_opt, V_W_opt, V_F_opt, V_G_opt = get_system_state(optimal_delta)
        P_W_opt_pu = S_W_opt.real
        P_W_opt_mw = P_W_opt_pu * S_base
        
        final_pf = P_W_opt_pu / np.abs(S_W_opt)
        
        # Calculate system losses
        I_WF_opt = (V_W_opt - V_F_opt) * Y_WF
        I_GF_opt = (V_F_opt - V_G_opt) * Y_GF
        loss_WF_pu = np.abs(I_WF_opt)**2 * Z_WF.real
        loss_GF_pu = np.abs(I_GF_opt)**2 * Z_GF.real
        loss_F_pu = np.abs(V_F_opt)**2 / R_F
        fundamental_loss_pu = loss_WF_pu + loss_GF_pu + loss_F_pu
        total_loss_pu = fundamental_loss_pu * 1.06 # Add 6% for harmonics
        total_loss_mw = total_loss_pu * S_base
        
        print("Optimization Successful.\n")
        print("--- Optimal System State ---")
        print(f"Optimal Reactive Power Injection (Q_opt): {Q_opt_mvar:.4f} MVAR")
        print(f"Required Active Power Generation (P_W): {P_W_opt_mw:.4f} MW")
        print(f"Final Power Factor at Bus-W: {final_pf:.4f} lagging")
        print(f"Total System Losses (incl. harmonics): {total_loss_mw:.4f} MW\n")
        
        print("--- Final Power Injection Equation ---")
        print("The required complex power injection at Bus-W is S = P + jQ.")
        print("The values making up this equation are:")
        print(f"P = {P_W_opt_mw:.4f} MW  ({P_W_opt_pu:.4f} p.u.)")
        print(f"Q = {Q_opt_mvar:.4f} MVAR ({Q_opt_pu:.4f} p.u.)")
        
        return Q_opt_mvar

    else:
        print("Optimization failed.")
        print(solution.message)
        return None

if __name__ == '__main__':
    Q_optimal = solve_fault_compensation()
