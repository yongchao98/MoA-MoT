import math

def print_magnetic_field_expression():
    """
    This function prints the derived symbolic expression for the magnetic field H_z(x, z).
    """
    
    # Define the symbolic variables in the equation as strings
    Ha = "Ha"  # Applied magnetic field
    H0 = "H0"  # Penetration field H0 = Jc*d/pi
    w = "w"    # Half-width of the strips
    D = "D"    # Stacking interval
    x = "x"    # x-coordinate
    z = "z"    # z-coordinate
    pi = "pi"  # The mathematical constant pi
    
    # Construct the terms of the equation
    # Term 1: Coefficient part of the induced field
    # pi**2 * H0 * w**2 / D**2
    coeff_part = f"({pi}**2 * {H0} * {w}**2 / {D}**2)"
    
    # Term 2: Numerator of the fractional part
    # 1 - cosh(2*pi*x/D)*cos(2*pi*z/D)
    numerator_arg = f"(2 * {pi} * {x} / {D})"
    numerator_part = f"(1 - cosh({numerator_arg}) * cos(2 * {pi} * {z} / {D}))"

    # Term 3: Denominator of the fractional part
    # (cosh(2*pi*x/D) - cos(2*pi*z/D))**2
    denominator_part = f"(cosh({numerator_arg}) - cos(2 * {pi} * {z} / {D}))**2"
    
    # Assemble the full expression for H_z(x, z)
    expression = f"H_z({x}, {z}) = {Ha} + {coeff_part} * ({numerator_part} / {denominator_part})"
    
    print("The expression for the magnetic field H_z(x, z) is:")
    print(expression)
    
    # Extract the final answer content for the specified format
    final_answer_content = f"{Ha} + ({pi}**2 * {H0} * {w}**2 / {D}**2) * (1 - cosh(2*{pi}*{x}/{D})*cos(2*{pi}*{z}/{D})) / (cosh(2*{pi}*{x}/{D}) - cos(2*{pi}*{z}/{D}))**2"
    
    # The problem asks to return the answer in a specific format at the end.
    # While the code prints the readable version, we format the core mathematical expression for the final output.
    print("\nFormatted answer for submission:")
    print(f"<<<{final_answer_content}>>>")

if __name__ == "__main__":
    print_magnetic_field_expression()
