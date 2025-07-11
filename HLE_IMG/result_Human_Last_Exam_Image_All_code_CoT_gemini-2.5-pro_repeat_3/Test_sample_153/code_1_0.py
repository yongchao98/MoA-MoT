def solve_power_system_equations():
    """
    This function formulates and prints the equations for total real power delivered and voltage drop
    based on the provided power system parameters.
    """

    # --- Part 1: Total Real Power Delivered ---

    # Transformer efficiency
    eta_transformer = 0.98

    # Transmission line resistive losses lead to an efficiency of (1 - loss)
    line_loss_percent = 2.0
    eta_line = 1 - (line_loss_percent / 100.0)

    # The total efficiency is the product of individual efficiencies
    # eta_total = eta_transformer * eta_line
    eta_total = eta_transformer * eta_line

    # The equation for delivered power is P_delivered = eta_total * (P_wind + P_pv)
    # We print the symbolic equation with the calculated coefficient.
    p_delivered_coeff = eta_total
    
    print("Equation for Total Real Power Delivered:")
    # We use .4f to ensure the output is 0.9604
    print(f"P_delivered = {p_delivered_coeff:.4f} * (P_wind + P_pv)")
    print("-" * 30)

    # --- Part 2: Voltage Drop Equation ---

    # The problem asks to identify the correct formula from the choices.
    # Based on our analysis, Option C is the most plausible.
    # We will reconstruct the voltage drop equation from Option C.
    
    # Real power coefficient (should match p_delivered_coeff)
    p_term_coeff = 0.9604

    # Reactive power compensation term coefficient from Option C
    q_comp_coeff = 0.979

    # Line impedance Z_line = R + jX
    R_line = 0.08
    X_line = 0.16

    print("Equation for Voltage Drop:")
    # We print the symbolic equation with all its numeric parts.
    print(f"V_drop = (({p_term_coeff:.4f} * (P_wind + P_pv)) + j * Q_comp * {q_comp_coeff}) / V_nominal * ({R_line} + j{X_line})")

solve_power_system_equations()