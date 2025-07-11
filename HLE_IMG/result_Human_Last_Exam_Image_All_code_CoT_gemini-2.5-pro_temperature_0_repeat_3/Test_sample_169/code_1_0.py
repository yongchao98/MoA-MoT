import math

def solve_wind_compensation():
    """
    Calculates the optimal reactive power injection to restore voltage after a fault.
    """
    # --- 1. Define System Parameters ---
    V_W_target_kV = 0.575  # Target voltage at Bus-W in kV
    R_WF_ohm = 0.01        # Resistance from Bus-W to Fault in Ohms
    X_WF_ohm = 0.05        # Reactance from Bus-W to Fault in Ohms
    R_F_ohm = 0.1          # Fault resistance in Ohms
    PF_min = 0.95          # Minimum lagging power factor of the generator
    Q_comp_max_MVAR = 10   # Maximum reactive power of the compensator in MVAR
    
    # Harmonic losses: 2% (3rd) + 4% (5th) = 6% increase in total losses
    harmonic_loss_increase = 0.06
    loss_factor = 1 + harmonic_loss_increase

    print("--- System Parameters ---")
    print(f"Target Voltage at Bus-W: {V_W_target_kV} kV")
    print(f"Transmission Impedance (Z_WF): {R_WF_ohm} + j{X_WF_ohm} Ohm")
    print(f"Fault Resistance (R_F): {R_F_ohm} Ohm")
    print(f"Generator Minimum Power Factor: {PF_min} lagging")
    print(f"Harmonic Loss Increase: {harmonic_loss_increase*100}%\n")

    # --- 2. Simplify the Model ---
    # The total impedance of the simplified load (line + fault)
    R_total_ohm = R_WF_ohm + R_F_ohm
    X_total_ohm = X_WF_ohm
    print("--- Simplified Model Calculations ---")
    print(f"Total series resistance (R_total): {R_WF_ohm} + {R_F_ohm} = {R_total_ohm:.3f} Ohm")
    print(f"Total series reactance (X_total): {X_total_ohm:.3f} Ohm\n")

    # --- 3. Calculate Required Power at Bus-W ---
    # Convert target voltage to Volts
    V_W_target_V = V_W_target_kV * 1000
    
    # Calculate the squared magnitude of the total impedance
    Z_total_sq_ohm2 = R_total_ohm**2 + X_total_ohm**2
    
    # Calculate the required active power (P_W) at Bus-W, including harmonic losses.
    # The formula is derived from P_W = loss_factor * (|I|^2 * R_total) and |S|^2 = |V|^2 * |I|^2
    # which leads to P_W = (V_W^2 * R_total) / (loss_factor * Z_total^2)
    P_W_W = (V_W_target_V**2 * R_total_ohm) / (loss_factor * Z_total_sq_ohm2)
    P_W_MW = P_W_W / 1e6

    # The total reactive power (Q_W) is determined by the impedance ratio
    Q_W_MVAR = P_W_MW * (X_total_ohm / R_total_ohm)
    
    print("--- Power Requirements at Bus-W (to maintain target voltage) ---")
    print(f"Required Active Power (P_W): {P_W_MW:.3f} MW")
    print(f"Required Reactive Power (Q_W): {Q_W_MVAR:.3f} MVAR\n")

    # --- 4. Apply Generator Constraints ---
    # The generator must supply all the active power
    P_gen_MW = P_W_MW
    
    # The maximum reactive power the generator can supply is limited by its power factor
    # Q_gen_max = P_gen * tan(phi), where cos(phi) = PF_min
    tan_phi_max = math.tan(math.acos(PF_min))
    Q_gen_max_MVAR = P_gen_MW * tan_phi_max
    
    print("--- Generator Contribution ---")
    print(f"Active Power from Generator (P_gen): {P_gen_MW:.3f} MW")
    print(f"Max Reactive Power from Generator (Q_gen_max) at PF={PF_min}: {Q_gen_max_MVAR:.3f} MVAR\n")

    # --- 5. Solve for Optimal Compensator Injection ---
    # The compensator must supply the remaining reactive power.
    # To minimize Q_comp, the generator should supply its maximum possible Q.
    Q_opt_MVAR = Q_W_MVAR - Q_gen_max_MVAR
    
    print("--- Optimal Reactive Power Injection (Q_opt) ---")
    print("The optimal (minimum) reactive power from the compensator is the total required reactive power minus the maximum reactive power the generator can provide.")
    print(f"Q_opt (MVAR) = Q_W (MVAR) - Q_gen_max (MVAR)")
    print(f"Q_opt (MVAR) = {Q_W_MVAR:.3f} - {Q_gen_max_MVAR:.3f}")
    
    # Final result
    print(f"\nFinal Answer: The optimal reactive power injection required is {Q_opt_MVAR:.3f} MVAR.")

    # Check if the solution is feasible
    if Q_opt_MVAR > Q_comp_max_MVAR:
        print(f"Warning: Required injection ({Q_opt_MVAR:.3f} MVAR) exceeds compensator capacity ({Q_comp_max_MVAR} MVAR).")
    elif Q_opt_MVAR < 0:
        print("The generator can supply all required reactive power. No compensation is needed (Q_opt = 0 MVAR).")
    else:
        print(f"This is within the compensator's capacity of {Q_comp_max_MVAR} MVAR.")
        
    return Q_opt_MVAR

# Run the solver
optimal_q = solve_wind_compensation()
# The final answer format requested by the user prompt
# print(f"\n<<<{optimal_q:.3f}>>>")