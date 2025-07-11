import math

def print_magnetic_field_expression():
    """
    This function prints the derived expressions for the magnetic field components
    Hz and Hx for a stack of superconducting strips under the condition |x| >> a.
    """

    # --- Introduction ---
    print("This script provides the expression for the magnetic field (H) for a stack")
    print("of thin superconducting strips in an applied field Ha, under the condition |x| >> a.")
    print("The solution is based on the model by Brandt and Zeldov.")
    print("-" * 70)

    # --- Definitions ---
    print("Variable Definitions:")
    print("  H_z, H_x: Components of the magnetic field in the z and x directions.")
    print("  x, z: Spatial coordinates.")
    print("  Jc: Critical current density of the superconductor.")
    print("  d: Thickness of each strip.")
    print("  D: Stacking interval (distance between the centers of strips).")
    print("  a: Position of the flux front, determined by the applied field Ha.")
    print("  pi: The mathematical constant pi.")
    print("-" * 70)

    # --- Characteristic Field H0 ---
    print("The characteristic field scale, H0, is defined as:")
    print("  H0 = (Jc * d) / pi")
    print("-" * 70)
    
    # --- Final Expressions ---
    print("For |x| >> a, the magnetic field components are approximated by:")

    # Expression for Hz(x, z)
    hz_expression = "H0 * (ln(2) - ln(sinh(pi*a/D)) + 0.5*ln(sinh(pi*x/D)**2 + sin(pi*z/D)**2))"
    print("\nThe z-component of the magnetic field is:")
    print(f"  Hz(x, z) ≈ {hz_expression}")

    # Expression for Hx(x, z)
    # The atan2 function is used for numerical stability and correctness across all quadrants.
    # atan2(y, x) computes atan(y/x)
    hx_expression = "H0 * atan2(tan(pi*z/D), tanh(pi*x/D))"
    print("\nThe x-component of the magnetic field is:")
    print(f"  Hx(x, z) ≈ {hx_expression}")

    print("-" * 70)
    print("Note: 'ln' is the natural logarithm, and 'atan2(y,x)' is the arc tangent of y/x.")

# Execute the function to print the result
print_magnetic_field_expression()