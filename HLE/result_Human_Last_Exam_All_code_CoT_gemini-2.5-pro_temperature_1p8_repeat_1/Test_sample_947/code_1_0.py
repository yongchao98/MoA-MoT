import math

def print_magnetic_field_formula():
    """
    Prints the symbolic formula for the magnetic field in a stack of superconducting strips.
    """
    
    # Define variables as strings for symbolic representation
    Jc = "Jc"  # Critical current density
    d = "d"    # Thickness of the strip
    pi = "π"   # The mathematical constant pi
    H0 = "H0"  # Characteristic field from the problem
    x = "x"    # Position along the strip width
    w = "w"    # Half-width of the strip
    D = "D"    # Stacking interval

    # The expression for the magnetic field H_z(x) in the plane of the strips (z=0)
    # is derived from the fully penetrated critical state model for an infinite stack.
    # The condition |x| >> a allows us to approximate the penetration depth a -> 0.
    
    # Construct the terms of the formula as strings
    term1 = f"tan({pi}*{x}/{D})"
    term2 = f"tan({pi}*{w}/{D})"
    
    numerator = f"({term1} - {term2})"
    denominator = f"({term1} + {term2})"
    
    # The full formula string
    # Using H0 as defined in the problem description
    final_formula = f"H_z({x}) = {H0} * ln| {numerator} / {denominator} |"
    
    # Print the final result
    print("Based on the fully penetrated critical state model for an infinite stack of strips, the expression for the magnetic field is:")
    print(final_formula)
    print("\nWhere the variables represent:")
    print(f"  H_z({x}): Magnetic field in the z-direction at position x.")
    print(f"  {H0}: Characteristic field, defined as {Jc} * {d} / {pi}.")
    print(f"  {Jc}: Critical current density of the superconductor.")
    print(f"  {d}: Thickness of each strip.")
    print(f"  {w}: Half-width of each strip (total width is 2w).")
    print(f"  {D}: The stacking period or interval between strips.")
    print(f"  {x}: The position coordinate in the direction of the width.")
    print(f"  {pi}: The mathematical constant pi (approx. 3.14159).")
    print( "  ln|...|: The natural logarithm of the absolute value of the argument.")

# Execute the function to print the formula
print_magnetic_field_formula()
