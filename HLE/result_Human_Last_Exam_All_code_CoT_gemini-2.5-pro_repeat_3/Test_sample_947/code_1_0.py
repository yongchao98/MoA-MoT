import math

def print_magnetic_field_expression():
    """
    This function prints the derived mathematical expressions for the magnetic field
    components (Hx and Hz) for a stack of superconducting strips.
    """

    # --- Explanation of variables ---
    # H_a: Applied magnetic field in the z-direction
    # H_0: Full penetration field, given by Jc * d / pi
    # w: Half-width of each superconducting strip
    # D: Stacking interval (distance between the centers of adjacent strips)
    # x: Coordinate in the direction of the strip's width
    # z: Coordinate in the direction of stacking
    # pi: Mathematical constant pi

    # Using unicode for pi for better readability
    pi = "\u03C0"

    # --- Constructing the expression strings ---

    # Common pre-factor in the field expressions from the strips
    # This term arises from the magnetic moment and the geometry of the array.
    # It is (pi^2 * H0 * w^2) / (2 * D^2)
    pre_factor = f"({pi}**2 * H_0 * w**2) / (2 * D**2)"

    # Common denominator in both field components
    # It is (cosh(2*pi*x/D) - cos(2*pi*z/D))^2
    denominator = f"(cosh(2*{pi}*x/D) - cos(2*{pi}*z/D))**2"

    # Numerator for the H_x component
    # It is 1 - cosh(2*pi*x/D) * cos(2*pi*z/D)
    hx_numerator = f"(1 - cosh(2*{pi}*x/D) * cos(2*{pi}*z/D))"

    # Numerator for the H_z component from the strips
    # It is sinh(2*pi*x/D) * sin(2*pi*z/D)
    hz_stack_numerator = f"sinh(2*{pi}*x/D) * sin(2*{pi}*z/D)"

    # --- Final Expressions ---

    # Expression for the x-component of the magnetic field
    hx_expression = f"H_x(x, z) = {pre_factor} * {hx_numerator} / {denominator}"

    # Expression for the z-component of the magnetic field
    # It's the sum of the applied field and the field from the strips
    hz_expression = f"H_z(x, z) = H_a - {pre_factor} * {hz_stack_numerator} / {denominator}"

    print("The expression for the x-component of the magnetic field is:")
    print(hx_expression)
    print("\nThe expression for the z-component of the magnetic field is:")
    print(hz_expression)

# Execute the function to print the expressions
print_magnetic_field_expression()