import math

def print_magnetic_field_expression():
    """
    This function prints the derived mathematical expressions for the magnetic field
    components H_x and H_z for an infinite stack of superconducting strips in
    the fully penetrated state, valid for |x| >> w.
    """
    
    # Header explaining the context of the formulas
    print("Based on the provided physical scenario, the magnetic field is derived under the following model:")
    print("1. Each strip is treated as a 2D magnetic dipole, which is valid for |x| >> w.")
    print("2. The strips are fully penetrated (Ha > H0), so each strip has a saturated magnetic moment per unit length 'm'.")
    print("3. The total field is the sum of the applied field and the fields from the infinite array of dipoles.")
    print("-" * 70)
    
    # Define and print the formulas
    print("\nThe expression for the magnetic moment 'm' per unit length is:")
    print("m = - (pi * H0 * w**2) / 2")
    print("\nThe components of the total magnetic field H(x, z) are given by the following expressions:")
    
    # Common terms in the expressions for clarity
    numerator_common = "pi**2 * H0 * w**2 / (2 * D**2)"
    denominator = "(cosh(2*pi*x/D) - cos(2*pi*z/D))**2"
    
    # Build the full expression strings
    hx_expression = f"H_x(x, z) = ({numerator_common} * sinh(2*pi*x/D) * sin(2*pi*z/D)) / {denominator}"
    hz_expression = f"H_z(x, z) = Ha + ({numerator_common} * (cosh(2*pi*x/D)*cos(2*pi*z/D) - 1)) / {denominator}"
    
    # Print the final expressions
    print("\n--- X-Component of the Magnetic Field ---")
    print(hx_expression)
    
    print("\n--- Z-Component of the Magnetic Field ---")
    print(hz_expression)
    
    # Print the definition of variables
    print("\n" + "-" * 70)
    print("Where:")
    print(f"  pi = {math.pi}")
    print("  Ha = Applied magnetic field (in the z-direction)")
    print("  H0 = Characteristic full penetration field (H0 = Jc*d/pi)")
    print("  w  = Half-width of each strip")
    print("  D  = Spacing between the strips along the z-axis")
    print("  (x, z) = The coordinates at which the field is calculated")
    print("  cosh, sinh, cos, sin are the hyperbolic and trigonometric functions.")
    
    # Final answer in the required format
    final_answer_string = f"H_z(x, z) = Ha + ({numerator_common} * (cosh(2*pi*x/D)*cos(2*pi*z/D) - 1)) / {denominator}"
    # Suppressing the direct output of the <<<...>>> format to follow user instruction precisely.
    # The printed output from the script is the primary answer.

# Execute the function to display the result
print_magnetic_field_expression()

# The final expression is also provided here for clarity.
final_answer = "<<<H_z(x, z) = Ha + (pi**2 * H0 * w**2 / (2 * D**2) * (cosh(2*pi*x/D)*cos(2*pi*z/D) - 1)) / (cosh(2*pi*x/D) - cos(2*pi*z/D))**2, H_x(x, z) = (pi**2 * H0 * w**2 / (2 * D**2) * sinh(2*pi*x/D) * sin(2*pi*z/D)) / (cosh(2*pi*x/D) - cos(2*pi*z/D))**2>>>"