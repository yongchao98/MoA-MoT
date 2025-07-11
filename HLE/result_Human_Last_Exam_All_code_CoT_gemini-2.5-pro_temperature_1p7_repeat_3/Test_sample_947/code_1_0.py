import math

def print_magnetic_field_formula():
    """
    This function prints the step-by-step derivation and the final formula
    for the magnetic field of a stack of superconducting strips.
    """
    
    print("The total magnetic field H_z(x, z) is the superposition of the applied field (Ha) and the field induced by the currents in the strips (H_strips).")
    print("The formula is derived assuming the strips are in the fully penetrated critical state, as implied by Ha > H0.")
    print("-" * 50)
    
    # --- Final Equation Structure ---
    print("Final Equation Structure:")
    print("H_z(x, z) = Ha + H_strips(x, z)")
    print("-" * 50)

    # --- H_strips Term ---
    print("The term for the field from the strips is:")
    print("H_strips(x, z) = (Jc * d / (4 * pi)) * ln(F(x, z))")
    print("\nWhere:")
    print("  Jc: Critical current density")
    print("  d: Thickness of each strip")
    print("  pi: The constant pi (approx 3.14159)")
    print("  ln: The natural logarithm")
    print("  F(x, z): A geometrical factor depending on position (x, z) and strip geometry (w, D)")
    print("-" * 50)
    
    # --- Geometrical Factor F(x, z) ---
    print("The geometrical factor F(x, z) is given by:")
    print("F(x, z) = ( A(x, z)**2 ) / ( A(x+w, z) * A(x-w, z) )")
    print("\nWhere:")
    print("  w: Half-width of each strip (total width is 2w)")
    print("  A(u, z): A function that arises from summing the fields of the infinite stack.")
    print("-" * 50)

    # --- Function A(u, z) ---
    print("The function A(u, z) is defined as:")
    print("A(u, z) = sinh(pi * u / D)**2 + sin(pi * z / D)**2")
    print("\nWhere:")
    print("  sinh: The hyperbolic sine function")
    print("  sin: The sine function")
    print("  D: The stacking interval between strips")
    print("  u: A placeholder for the horizontal position (e.g., x, x+w, x-w)")
    print("-" * 50)
    
    # --- Full Expression ---
    print("Combining all the parts, the complete expression for the total magnetic field is:")
    final_formula = "H_z(x, z) = Ha + (Jc*d / (4*pi)) * ln[ (sinh(pi*x/D)**2 + sin(pi*z/D)**2)**2 / ((sinh(pi*(x+w)/D)**2 + sin(pi*z/D)**2) * (sinh(pi*(x-w)/D)**2 + sin(pi*z/D)**2)) ]"
    print(final_formula)

if __name__ == '__main__':
    print_magnetic_field_formula()
<<<H_z(x, z) = Ha + (Jc*d / (4*pi)) * ln[ (sinh(pi*x/D)**2 + sin(pi*z/D)**2)**2 / ((sinh(pi*(x+w)/D)**2 + sin(pi*z/D)**2) * (sinh(pi*(x-w)/D)**2 + sin(pi*z/D)**2)) ]>>>