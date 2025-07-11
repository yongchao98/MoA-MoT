import numpy as np

def solve_power_injection():
    """
    Formulates and solves the problem to find the optimal reactive power injection.
    """
    # 1. Define system parameters in per-unit (p.u.)
    # Given transmission impedance from the problem description
    R_wf = 0.01  # p.u.
    X_wf = 0.05  # p.u.

    # Fault resistance is given as 0.1 Ohm. We must assume it is already in p.u.
    # as other interpretations lead to physically impossible scenarios.
    R_f = 0.1  # p.u.

    # Target voltage for Bus-W
    V_target = 0.575  # p.u.

    # System base power
    S_base = 100  # MVA
    
    # Maximum reactive power capacity of the compensating device
    Q_max = 10    # MVAR

    # Harmonic losses increase total system losses by 6%
    harmonic_loss_factor = 1.06

    print("--- Step 1: System Parameter Setup (in per-unit) ---")
    print(f"Base transmission impedance: Z_WF = {R_wf} + j{X_wf} p.u.")
    print(f"Fault resistance: R_f = {R_f} p.u.")
    print(f"Target voltage at Bus-W: V_target = {V_target} p.u.")
    print(f"Harmonic loss factor: {harmonic_loss_factor-1.0:.0%}\n")

    # 2. Adjust for Harmonic Losses
    # The 6% increase in losses is applied to the resistive part of the impedance.
    R_eff = R_wf * harmonic_loss_factor
    Z_sq = R_eff**2 + X_wf**2
    Z_abs = np.sqrt(Z_sq)

    print("--- Step 2: Incorporate Harmonic Losses ---")
    print(f"Effective line resistance, R_eff = {R_wf} * {harmonic_loss_factor} = {R_eff:.4f} p.u.")
    print(f"Effective line impedance, Z_eff = {R_eff:.4f} + j{X_wf:.4f} p.u. (Magnitude |Z| = {Z_abs:.4f} p.u.)\n")

    # 3. Solve for the required power angle (delta)
    # With the wind farm output Pw=0, real power from the grid must supply the fault.
    # P_in(from grid) = P_fault = V_target^2 / R_f
    # This sets up a non-linear equation for the angle delta.
    P_fault = V_target**2 / R_f
    
    # The equation is: X_wf * sin(delta) - R_eff * cos(delta) = C
    C = V_target * (Z_sq / R_f - R_eff)
    
    # Solve A*sin(d) + B*cos(d) = C_prime form using atan2 for stability
    # Our equation: X_wf*sin(d) - R_eff*cos(d) = C
    A_trig = X_wf
    B_trig = -R_eff
    
    alpha = np.arctan2(B_trig, A_trig)
    sin_delta_minus_alpha = C / Z_abs
    
    # Ensure the solution is physically possible (arcsin input is between -1 and 1)
    if abs(sin_delta_minus_alpha) > 1:
        print("Error: Physical solution does not exist. The line cannot supply the required fault power at this voltage.")
        return
        
    delta = np.arcsin(sin_delta_minus_alpha) + alpha

    print("--- Step 3: Solve for Power Angle (delta) ---")
    print("Equation from real power balance (P_in = P_fault):")
    final_eq_for_delta = f"{X_wf:.2f}*sin(delta) - {R_eff:.4f}*cos(delta) = {V_target:.3f} * ({Z_sq:.6f} / {R_f:.1f} - {R_eff:.4f})"
    print(final_eq_for_delta)
    print(f"Which simplifies to: {X_wf:.2f}*sin(delta) - {R_eff:.4f}*cos(delta) = {C:.5f}")
    print(f"Solved angle, delta = {np.rad2deg(delta):.2f} degrees ({delta:.4f} radians)\n")
    
    # 4. Calculate the reactive power injection Q_opt
    # From reactive power balance: Q_in + Q_opt = 0 => Q_opt = -Q_in
    # Using the standard power flow equation for Q_in
    q_in_bracket = V_target * X_wf - R_eff * np.sin(delta) - X_wf * np.cos(delta)
    Q_in = (V_target / Z_sq) * q_in_bracket
    Q_opt_pu = -Q_in
    Q_opt_mvar = Q_opt_pu * S_base
    
    print("--- Step 4: Calculate Optimal Reactive Power (Q_opt) ---")
    print("From reactive power balance, Q_opt = -Q_in. Q_in is calculated as:")
    q_in_eq_str_1 = f"Q_in = ( V_target / |Z|^2 ) * ( V_target*X_wf - R_eff*sin(delta) - X_wf*cos(delta) )"
    print(q_in_eq_str_1)
    
    print("Plugging in the numbers:")
    final_calc_eq = f"Q_opt = - ( {V_target:.3f} / {Z_sq:.6f} ) * ( {V_target*X_wf:.4f} - {R_eff*np.sin(delta):.4f} - {X_wf*np.cos(delta):.4f} )"
    print(final_calc_eq)
    final_calc_eq_simpl = f"Q_opt = - ( {V_target/Z_sq:.2f} ) * ( {q_in_bracket:.4f} )"
    print(final_calc_eq_simpl)
    
    print(f"Final calculated Q_opt = {Q_opt_pu:.4f} p.u.\n")

    print("--- Final Answer ---")
    print(f"The optimal reactive power injection Q_opt is {Q_opt_pu:.2f} p.u.")
    print(f"In physical units, this is {Q_opt_pu:.2f} * {S_base} MVA = {Q_opt_mvar:.2f} MVAR.")

solve_power_injection()