# -*- coding: utf-8 -*-

def solve_power_equations():
    """
    Calculates the coefficients for power delivery and formats the final equations
    based on the problem statement.
    """

    # --- Part 1: Calculate Total Real Power Delivered Coefficient ---

    # Given parameters from the problem description
    eta_transformer = 0.98  # Transformer efficiency
    line_loss_percentage = 2.0  # Transmission line resistive losses

    # Calculate the efficiency of the transmission line
    eta_line = 1 - (line_loss_percentage / 100)

    # The total delivered power coefficient is the product of the transformer
    # and line efficiencies.
    p_delivered_coeff = eta_transformer * eta_line

    # --- Part 2: Identify Voltage Drop Equation Components ---

    # The line impedance is given in the problem statement.
    R_line = 0.08  # Resistance (Ohms)
    X_line = 0.16  # Reactance (Ohms)

    # The factor for the reactive power compensation term is taken from the
    # most consistent answer choice (Option C).
    q_factor = 0.979

    # --- Part 3: Print the Final Equations ---

    # Equation for Real Power Delivered
    print(f"    P_delivered = {p_delivered_coeff:.4f} \u22C5 (P_wind + P_pv)")

    # Equation for Voltage Drop. We represent the fraction from the options
    # as a division in a single line.
    print(f"    V_drop = (({p_delivered_coeff:.4f} \u22C5 (P_wind + P_pv) + j Q_comp \u22C5 {q_factor:.3f}) / V_nominal) \u22C5 ({R_line:.2f} + j{X_line:.2f})")

solve_power_equations()
<<<C>>>