import math

def display_magnetization_formula():
    """
    This function prints the analytical expression for the initial magnetization
    curve of a superconducting strip based on the critical-state model.
    """

    # --- Variable Definitions ---
    # We define the symbols used in the equation for clarity.
    M_H = "M(H)"  # Magnetization as a function of H
    H = "H"       # Applied magnetic field
    Jc = "Jc"     # Critical current density
    b = "b"       # Half-thickness of the superconducting strip

    # The full penetration field, Hp, is the field required to penetrate
    # to the center of the slab.
    Hp_formula = f"{Jc} * {b}"
    Hp_symbol = "Hp"

    # --- Initial Magnetization Curve Expression ---
    print("Based on the critical-state model for a wide, thin superconducting strip in a parallel applied magnetic field, the analytical expression for the initial magnetization curve is derived.")
    print("\n--- Key Variables ---")
    print(f"M(H): Magnetization")
    print(f"H: Applied magnetic field")
    print(f"Jc: Critical current density")
    print(f"b: Half-thickness of the strip")
    print(f"Hp: Full penetration field, where Hp = {Hp_formula}")

    print("\n--- Initial Magnetization Curve M(H) ---")

    # Case 1: Partial Field Penetration
    print(f"\n1. For partial penetration (0 <= {H} <= {Hp_symbol}):")
    # In this regime, the magnetization is a quadratic function of the applied field.
    # The equation is M(H) = H^2 / (2 * b * Jc) - H
    numerator = f"{H}^2"
    denominator = f"2 * {b} * {Jc}"
    equation_part1 = f"{M_H} = ({numerator}) / ({denominator}) - {H}"
    print(f"   {equation_part1}")
    print("   Each variable in the equation is printed above.")


    # Case 2: Full Field Penetration
    print(f"\n2. For full penetration ({H} > {Hp_symbol}):")
    # In this regime, the superconductor is saturated, and the magnetization is constant.
    # The equation is M(H) = - (Jc * b) / 2
    numerator_2 = f"{Jc} * {b}"
    denominator_2 = "2"
    equation_part2 = f"{M_H} = -({numerator_2}) / {denominator_2}"
    print(f"   {equation_part2}")
    print("   Each variable in the equation is printed above.")


if __name__ == '__main__':
    display_magnetization_formula()
