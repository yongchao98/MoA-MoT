import math

def solve_magnetization():
    """
    Calculates and prints the analytical expression for the initial magnetization
    curve of a superconducting slab based on the critical-state model.
    """

    # --- Problem Parameters ---
    # These are example values. The final expression is analytical.
    # 'a' is the half-thickness of the slab in meters.
    a = 0.001  # e.g., 1 mm
    # 'Jc' is the critical current density in Amperes per square meter (A/m^2).
    Jc = 1e7   # e.g., 10^7 A/m^2

    # --- Derived Quantity ---
    # H_p is the full penetration field, the field at which the shielding
    # currents reach the center of the slab.
    H_p = a * Jc

    # --- Output the Explanation and Solution ---
    print("This solution models the superconducting bar as an infinite slab of half-thickness 'a'.")
    print("The applied magnetic field 'H' is parallel to the slab surface.")
    print("The critical-state model with a constant critical current 'Jc' is used.")
    print("\n--- Analytical Expressions for Initial Magnetization M(H) ---\n")
    print(f"Given parameters:")
    print(f"  a (half-thickness) = {a} m")
    print(f"  Jc (critical current density) = {Jc:.0e} A/m^2\n")

    # The full penetration field, H_p = a * Jc
    print("The characteristic full penetration field H_p is calculated as:")
    print(f"  H_p = a * Jc = {a} * {Jc:.0e} = {H_p:.0f} A/m\n")

    # --- Regime 1: Partial Penetration (0 <= H <= H_p) ---
    print(f"1. For an applied field H in the range 0 <= H <= {H_p:.0f} A/m (Partial Penetration):")
    print("   The magnetization M is given by the formula:")
    # We output the literal formula with each component shown
    print(f"   M(H) = -H * (1 - H / (2 * a * Jc))")
    print(f"   Substituting the values for a and Jc:")
    print(f"   M(H) = -H * (1 - H / (2 * {a} * {Jc:.0e}))")
    print(f"   M(H) = -H * (1 - H / {2 * a * Jc}) \n")


    # --- Regime 2: Full Penetration (H > H_p) ---
    print(f"2. For an applied field H > {H_p:.0f} A/m (Full Penetration):")
    print("   The magnetization M saturates to a constant value:")
    # We output the literal formula with each component shown
    print(f"   M(H) = - (a * Jc) / 2")
    print(f"   Or simply M(H) = -H_p / 2")
    # Calculating the saturated value
    M_sat = -H_p / 2
    print(f"   Substituting the values:")
    print(f"   M(H) = -({a} * {Jc:.0e}) / 2 = {M_sat:.0f} A/m")


# Execute the function to print the solution
solve_magnetization()

# The final answer is the two-part expression for M(H)
final_expression_part1 = "For 0 <= H <= a*Jc: M(H) = -H * (1 - H / (2 * a * Jc))"
final_expression_part2 = "For H > a*Jc: M(H) = -a*Jc / 2"
# Combining them for the final answer format
final_answer = f"({final_expression_part1}); ({final_expression_part2})"
print(f"\n<<<{final_answer}>>>")