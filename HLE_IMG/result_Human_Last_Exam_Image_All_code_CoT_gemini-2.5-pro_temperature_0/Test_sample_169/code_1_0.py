import numpy as np
from scipy.optimize import minimize

def solve_reactive_power_optimization():
    """
    Calculates the optimal reactive power injection to mitigate a fault
    in a wind energy system by solving a nonlinear optimization problem.
    """
    # --- 1. Define System Parameters ---
    # Given values from the problem description
    V_nominal_kV = 0.575  # Nominal voltage at Bus-W, kV
    S_base_MVA = 100.0    # Base power, MVA
    Z_WF_ohm = 0.01 + 0.05j # Transmission impedance from Bus-W to Fault, Ohm
    R_F_ohm = 0.1         # Fault resistance, Ohm
    Q_max_MVAR = 10.0     # Max reactive power of compensator, MVAR
    PF_min = 0.95         # Minimum lagging power factor for the source
    harmonic_loss_increase = 0.06 # 6% increase in losses due to harmonics

    # The target voltage for Bus-W is its nominal voltage
    V_target_kV = 0.575

    # --- 2. Per-Unit System Conversion ---
    # Calculate base values for the per-unit system
    V_base_V = V_nominal_kV * 1000
    S_base_VA = S_base_MVA * 1e6
    Z_base_ohm = V_base_V**2 / S_base_VA

    # Convert given physical values to per-unit (p.u.)
    V_W_pu = V_target_kV / V_nominal_kV # Target voltage is 1.0 p.u.
    Z_WF_pu = Z_WF_ohm / Z_base_ohm
    R_F_pu = R_F_ohm / Z_base_ohm
    Q_max_pu = Q_max_MVAR / S_base_MVA

    # --- 3. Power Calculation for the Fault Path ---
    # The total impedance from Bus-W to ground through the fault
    Z_total_pu = Z_WF_pu + R_F_pu

    # Apparent power drawn by the fault path to maintain V_W at 1.0 p.u.
    # S = V*I* = V * (V/Z)* = |V|^2 / Z*
    S_fault_pu = V_W_pu**2 / np.conj(Z_total_pu)

    P_fault_pu = S_fault_pu.real
    Q_fault_pu = S_fault_pu.imag

    # --- 4. Incorporate Harmonic Losses ---
    # Harmonic losses increase the total active power loss by the specified percentage
    P_harmonic_loss_pu = P_fault_pu * harmonic_loss_increase
    # Total active power required at Bus-W
    P_load_pu = P_fault_pu + P_harmonic_loss_pu
    # Total reactive power required at Bus-W is consumed by the line's reactance
    Q_load_pu = Q_fault_pu

    # --- 5. Formulate and Solve the Optimization Problem ---
    # The power from the source is S_source = P_load + j(Q_load - Q_comp)
    # We want to minimize Q_comp, which we define as x[0]
    objective_function = lambda x: x[0]

    # Constraint 1: Power Factor >= 0.95 lagging.
    # This can be expressed as: Q_source / P_source <= tan(arccos(PF_min))
    # (Q_load_pu - x[0]) / P_load_pu <= tan(arccos(PF_min))
    # Rearranging for the solver's required format (g(x) >= 0):
    # P_load_pu * tan(arccos(PF_min)) - (Q_load_pu - x[0]) >= 0
    tan_phi_max = np.tan(np.arccos(PF_min))
    cons1 = {'type': 'ineq', 'fun': lambda x: P_load_pu * tan_phi_max - (Q_load_pu - x[0])}

    # Constraint 2: Power Factor must be "lagging", meaning Q_source > 0.
    # Q_load_pu - x[0] >= 0
    cons2 = {'type': 'ineq', 'fun': lambda x: Q_load_pu - x[0]}

    # Combine all constraints
    constraints = [cons1, cons2]

    # Define the bounds for Q_comp_pu: 0 <= Q_comp <= Q_max
    bounds = [(0, Q_max_pu)]

    # Provide an initial guess for the solver
    x0 = [0.0]

    # Solve the optimization problem using Sequential Least Squares Programming (SLSQP)
    solution = minimize(objective_function, x0, method='SLSQP', bounds=bounds, constraints=constraints)

    Q_opt_pu = solution.x[0]
    Q_opt_MVAR = Q_opt_pu * S_base_MVA

    # --- 6. Print the Results Step-by-Step ---
    print("--- Step 1: System Parameters and Per-Unit Conversion ---")
    print(f"Base Power (S_base): {S_base_MVA} MVA")
    print(f"Base Voltage (V_base): {V_nominal_kV} kV")
    print(f"Base Impedance (Z_base): {Z_base_ohm:.6f} Ohms")
    print(f"Total Fault Path Impedance (Z_total): {Z_total_pu.real:.4f} + j{Z_total_pu.imag:.4f} p.u.")

    print("\n--- Step 2: Power Analysis at Bus-W ---")
    print(f"Power drawn by fault (S_fault): {P_fault_pu:.5f} + j{Q_fault_pu:.5f} p.u.")
    print(f"Harmonic Losses (P_harmonic): {P_harmonic_loss_pu:.5f} p.u. ({harmonic_loss_increase*100:.0f}%)")
    print(f"Total Load at Bus-W (P_load + jQ_load): {P_load_pu:.5f} + j{Q_load_pu:.5f} p.u.")

    print("\n--- Step 3: Optimization and Final Calculation ---")
    print("The optimization goal is to find the minimum Q_comp that satisfies the power factor constraint.")
    print(f"The constraint is: Q_comp >= Q_load - P_load * tan(arccos({PF_min}))")
    
    # Display the final equation with numbers
    print("\nFinal Equation:")
    print(f"Q_opt_p.u. = {Q_load_pu:.5f} - ({P_load_pu:.5f} * {tan_phi_max:.4f})")
    
    min_q_comp_for_pf = Q_load_pu - P_load_pu * tan_phi_max
    print(f"Q_opt_p.u. = {Q_load_pu:.5f} - {P_load_pu * tan_phi_max:.5f}")
    print(f"Q_opt_p.u. = {min_q_comp_for_pf:.5f}")

    print("\n--- Final Answer ---")
    print(f"The optimal reactive power injection required is:")
    print(f"Q_opt = {Q_opt_pu:.5f} p.u.")
    print(f"Q_opt = {Q_opt_MVAR:.4f} MVAR")
    
    # Final answer in the required format
    print(f"\n<<< {Q_opt_MVAR:.4f} >>>")

# Execute the function
solve_reactive_power_optimization()