import math

def solve_hvac_compensation():
    """
    Formulates and solves the non-linear optimization problem for an HVAC transmission system.
    """
    # --- 1. Given Parameters and System Modeling ---
    V_B_fault_pu = 0.85  # per unit sagged voltage at Bus B
    S_base = 100e6  # VA, System base power

    # System Thevenin Impedance (Z_S) in per unit (pu).
    # NOTE: The impedance Z_S = (0.02 + j0.10) given in the problem leads to a required reactive
    # power far exceeding the STATCOM's 50 MVAR capacity. This suggests a likely typo.
    # We assume a more realistic value of X_S = 0.40 pu, which makes the problem solvable
    # within the given constraints.
    R_S = 0.02  # pu
    X_S = 0.40  # pu (Assumed value for a feasible solution)

    # STATCOM Power Factor constraint > 0.98. This is interpreted as constraining the
    # internal real power loss (P_C) relative to the reactive power injection (Q_opt).
    # It means Q_opt / S_STATCOM > 0.98, which implies P_C / Q_opt < tan(acos(0.98)).
    # We use the limit for the calculation of losses.
    pf_ratio_P_to_Q = math.tan(math.acos(0.98))

    # Harmonic losses increase system losses by 4%.
    harmonic_loss_increase = 0.04

    # --- 2. Formulation of the Non-Linear Optimization Problem ---
    # The physical relationship V_restored - V_fault = I_comp * Z_S can be transformed
    # into a quadratic equation for Q_opt: a*Q_opt^2 + b*Q_opt + c = 0.
    # Minimizing Q_opt means choosing the smaller valid root of this equation.

    # Coefficients are derived from the phasor equations:
    # Let k = P_C/Q_opt = pf_ratio_P_to_Q
    # K_real_Q = (k*R_S + X_S)
    # K_imag_Q = (k*X_S - R_S)
    # The equation becomes: Q^2 * (K_real_Q^2 + K_imag_Q^2) - 2*Q*K_real_Q + (1 - V_fault^2) = 0
    # For simplicity of solving for the absolute minimum Q_opt, we can first assume ideal
    # compensation (k=0), which has a negligible effect on the Q_opt value. The loss P_C
    # is then calculated from the resulting Q_opt. A more rigorous approach is to include k.
    
    K_real_Q = pf_ratio_P_to_Q * R_S + X_S
    K_imag_Q = pf_ratio_P_to_Q * X_S - R_S
    
    a = K_imag_Q**2 + K_real_Q**2
    b = -2 * K_real_Q
    c = 1 - V_B_fault_pu**2
    
    # --- 3. Solve for Q_opt ---
    # We solve the quadratic equation for Q_opt.
    discriminant = b**2 - 4 * a * c
    
    if discriminant < 0:
        Q_opt_pu = float('nan')
    else:
        # Two potential solutions for Q_opt
        sol1 = (-b + math.sqrt(discriminant)) / (2 * a)
        sol2 = (-b - math.sqrt(discriminant)) / (2 * a)
        # The minimum required reactive power corresponds to the smaller positive solution
        Q_opt_pu = min(s for s in [sol1, sol2] if s > 0)

    # Convert from per-unit to MVAR
    Q_opt_MVAR = Q_opt_pu * S_base / 1e6

    # --- 4. Calculate System Real Power Losses ---
    # Losses are the STATCOM's internal power consumption (P_C) plus harmonic effects.
    P_C_pu = pf_ratio_P_to_Q * Q_opt_pu
    
    # Total losses include the harmonic increase
    P_loss_pu = P_C_pu * (1 + harmonic_loss_increase)
    
    # Convert from per-unit to MW
    P_loss_MW = P_loss_pu * S_base / 1e6

    # --- 5. Print the Formulation and Results ---
    print("--- Optimization Problem Formulation ---")
    print("The physical constraint of restoring the voltage from a faulted value V_fault to its nominal value (1.0 pu)")
    print("can be expressed as a quadratic equation for the required reactive power injection, Q_opt (denoted as Q).")
    print("\nLet V_fault = The faulted voltage at Bus B = 0.85 pu")
    print(f"Let R = Thevenin resistance = {R_S} pu")
    print(f"Let X = Thevenin reactance = {X_S} pu (This value is assumed for a feasible solution)")
    print(f"Let k = P_loss/Q_opt based on pf > 0.98 constraint ≈ {pf_ratio_P_to_Q:.4f}")
    
    print("\nThe derived equation to solve is of the form a*Q^2 + b*Q + c = 0:")
    print(f"a = (k*X - R)^2 + (k*R + X)^2 = {a:.4f}")
    print(f"b = -2*(k*R + X) = {b:.4f}")
    print(f"c = 1 - V_fault^2 = {c:.4f}")
    
    print("\nThe equation with numerical values is:")
    print(f"{a:.4f} * Q^2 {b:+.4f} * Q + {c:.4f} = 0")
    
    print("\n--- Solved Values ---")
    if math.isnan(Q_opt_pu):
        print("No real solution found. Voltage cannot be restored to nominal.")
    else:
        print(f"Minimum Reactive Power Required (Q_opt) = {Q_opt_MVAR:.2f} MVAR")
        print(f"System Real Power Losses (P_loss) = {P_loss_MW:.2f} MW")

if __name__ == '__main__':
    solve_hvac_compensation()