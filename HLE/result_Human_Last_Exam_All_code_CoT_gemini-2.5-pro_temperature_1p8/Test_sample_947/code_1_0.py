import math

def print_magnetic_field_formula():
    """
    This function prints the symbolic formula for the magnetic field H_z(x, z)
    for an infinite stack of superconducting strips in a perpendicular magnetic field.

    The formula is valid for an applied field H_a > H_0, where H_0 = J_c*d/pi,
    and for positions |x| >> a, where 'a' is the flux penetration distance.
    """

    # --- Define variables as strings for printing the formula ---
    H_a, J_c, d, w, D, x, z, a, pi, H_0 = \
        "H_a", "J_c", "d", "w", "D", "x", "z", "a", "pi", "H_0"

    # --- Print explanatory text and preliminary definitions ---
    print("This script provides the formula for the magnetic field H_z(x, z).")
    print("The system is an infinite stack of superconducting strips of width 2w, thickness d,")
    print("and critical current density J_c, separated by a distance D.")
    print("-" * 70)
    
    # --- Definition of H_0 ---
    print(f"The characteristic field, H_0, is defined as:")
    print(f"{H_0} = ({J_c} * {d}) / {pi}")
    print("\n")

    # --- Relation for penetration depth 'a' ---
    print(f"For an applied field {H_a} > {H_0}, the flux penetrates to a distance 'a' from the edge.")
    print("This penetration depth 'a' is determined by the following relation:")
    print(f"cosh({pi} * {a} / {D}) = cosh({pi} * {w} / {D}) / cosh({H_a} / {H_0})")
    print("-" * 70)

    # --- Final formula for H_z(x,z) ---
    print("The total magnetic field distribution H_z(x, z) is given by the expression:")

    # Define components of the formula for clarity
    numerator_str = f"sin(2*{pi}*{x}/{D}) * sinh(2*{pi}*{z}/{D})"
    denominator_a_str = f"cos(2*{pi}*{a}/{D}) - cos(2*{pi}*{x}/{D}) * cosh(2*{pi}*{z}/{D})"
    denominator_w_str = f"cos(2*{pi}*{w}/{D}) - cos(2*{pi}*{x}/{D}) * cosh(2*{pi}*{z}/{D})"

    print(f"\nH_z({x}, {z}) = {H_a} + ({J_c} * {d}) / (2 * {pi}) * [ term_1 - term_2 ]\n")
    print("Where:")
    print(f"  term_1 = atan( ({numerator_str}) / ({denominator_a_str}) )")
    print(f"  term_2 = atan( ({numerator_str}) / ({denominator_w_str}) )")
    
# Execute the function to print the formula
print_magnetic_field_formula()