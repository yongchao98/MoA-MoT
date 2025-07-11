import cmath
import math

def solve_wind_compensation():
    """
    Solves the nonlinear optimization problem for reactive power compensation.
    """
    # 1. System Parameters
    S_base = 100e6  # VA
    V_target_actual = 0.575e3 # V
    V_base = 0.575e3 # V, chosen based on target and compensator nominal voltage
    
    Z_WF_ohm = complex(0.01, 0.05) # Ohms
    R_F_ohm = 0.1 # Ohms
    Q_max_mvar = 10 # MVAR
    PF_target = 0.95 # lagging
    harmonic_loss_factor = 1.06 # 6% increase

    # 2. Convert to Per-Unit
    Z_base = V_base**2 / S_base
    V_W_pu = V_target_actual / V_base

    Z_WF_pu = Z_WF_ohm / Z_base
    R_F_pu = R_F_ohm / Z_base
    Q_max_pu = Q_max_mvar * 1e6 / S_base

    # 3. Formulate and solve the problem
    # Model: |S_W| = |V_W|^2 / |Z_total|
    # Z_total is the impedance from Bus-W to ground through the fault
    Z_total_pu = Z_WF_pu + R_F_pu
    
    # Calculate the magnitude of the required apparent power from Bus-W
    S_W_mag_pu = V_W_pu**2 / abs(Z_total_pu)

    # From the power factor constraint (PF = P/S), we have P = S * PF
    # P^2 + Q^2 = S^2 => (S*PF)^2 + Q^2 = S^2 => Q = S * sqrt(1 - PF^2)
    P_W_pu = S_W_mag_pu * PF_target
    Q_W_pu = S_W_mag_pu * math.sqrt(1 - PF_target**2)

    # 4. Determine Optimal Reactive Power (Q_opt)
    # Assuming the generator provides no reactive power (unity PF),
    # the compensator must provide all of Q_W.
    Q_opt_pu = Q_W_pu
    
    # Check if the required compensation is within the device's capacity
    is_feasible = Q_opt_pu <= Q_max_pu

    # 5. Calculate Losses (for completeness)
    # Fundamental losses in the line Z_WF
    P_loss_line_pu = Z_WF_pu.real * (S_W_mag_pu**2 / V_W_pu**2)
    
    # Total losses including harmonic effects
    Total_loss_pu = P_loss_line_pu * harmonic_loss_factor
    Total_loss_mw = Total_loss_pu * S_base / 1e6

    # 6. Output Results
    print("--- System Parameters (Per Unit) ---")
    print(f"Base Voltage: {V_base/1000:.3f} kV")
    print(f"Base Power: {S_base/1e6:.0f} MVA")
    print(f"Base Impedance: {Z_base:.6f} Ohms")
    print(f"Target Voltage (V_W): {V_W_pu:.3f} p.u.")
    print(f"Transmission Impedance (Z_WF): {Z_WF_pu:.3f} p.u.")
    print(f"Fault Resistance (R_F): {R_F_pu:.3f} p.u.")
    print(f"Max Compensator Q (Q_max): {Q_max_pu:.3f} p.u.")
    print("-" * 35)
    
    print("\n--- Solved Power Flow at Bus-W ---")
    print(f"Required Apparent Power Magnitude (|S_W|): {S_W_mag_pu:.4f} p.u.")
    print(f"Required Active Power (P_W): {P_W_pu:.4f} p.u.")
    print(f"Required Reactive Power (Q_W): {Q_W_pu:.4f} p.u.")

    print("\n--- Final Equation with Solved Values ---")
    print("P_W^2 + Q_W^2 = |S_W|^2")
    print(f"{P_W_pu:.4f}^2 + {Q_W_pu:.4f}^2 = {S_W_mag_pu:.4f}^2")
    print(f"{P_W_pu**2:.4f} + {Q_W_pu**2:.4f} = {S_W_mag_pu**2:.4f}")
    
    print("\n--- Optimal Reactive Power Injection ---")
    if is_feasible:
        print(f"Optimal Reactive Power (Q_opt): {Q_opt_pu:.4f} p.u.")
        print(f"This is equivalent to {Q_opt_pu * S_base / 1e6:.4f} MVAR.")
        print(f"This is within the device's capacity of {Q_max_mvar} MVAR.")
    else:
        print(f"Required reactive power ({Q_opt_pu * S_base / 1e6:.2f} MVAR) exceeds device capacity ({Q_max_mvar} MVAR).")
        print("The voltage cannot be restored to the target under the given constraints.")

    print("\n--- System Losses ---")
    print(f"Total system losses (including harmonic effects): {Total_loss_mw:.4f} MW")

solve_wind_compensation()