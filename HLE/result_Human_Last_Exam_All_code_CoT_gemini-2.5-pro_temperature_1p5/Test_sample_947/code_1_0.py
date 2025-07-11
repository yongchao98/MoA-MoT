import math

def print_magnetic_field_expression():
    """
    This function prints the derived expression for the magnetic field H_z(x, z)
    for a stack of superconducting strips.
    The expression is valid for |x| >> a.
    """

    print("The expression for the z-component of the magnetic field H_z(x, z) is composed of the applied field H_a and a correction term from the induced currents in the strips.")
    print("\nThe general form is:")
    print("H_z(x, z) = H_a - C_amplitude * F_spatial(x, z)")
    print("\n-------------------------------------------------------------")

    # Define and print each component of the equation
    # This fulfills the "output each number in the final equation" requirement
    # by clearly showing all coefficients and variables.
    print("The complete expression for the magnetic field is:")
    print("H_z(x, z) = H_a - [ (π² * H₀ * w²) / (2 * D²) ] * tanh²(H_a / H₀) * [ cos(2πz / D) / sinh²(π|x| / D) ]")
    print("\nBreaking this down into its constituent parts:")
    
    # Applied Field
    print("\n1. Applied Field:")
    print("   H_a")

    # Amplitude part
    print("\n2. Amplitude Coefficient (C_amplitude):")
    print("   (π² * H₀ * w²) / (2 * D²) * tanh²(H_a / H₀)")
    print("   This term combines the fundamental constants and material properties (H₀), the geometry (w, D), and the dependence on the applied field (H_a).")
    
    # Spatial part
    print("\n3. Spatial Dependence (F_spatial):")
    print("   cos(2πz / D) / sinh²(π|x| / D)")
    print("   This term describes how the field correction varies in space. It is periodic along the stacking direction (z) and decays exponentially away from the stack in the x-direction.")

    # Variable definitions
    print("\n-------------------------------------------------------------")
    print("Variable definitions:")
    print("  H_z: z-component of the magnetic field at position (x, z)")
    print("  H_a: Applied external magnetic field")
    print("  H₀: Characteristic field of the superconductor, H₀ = Jc * d / π")
    print("  w: Half-width of the superconducting strips")
    print("  D: Stacking interval (distance between strips)")
    print("  Jc: Critical current density")
    print("  d: Thickness of the strips")
    print("  π: Mathematical constant pi (approx 3.14159)")
    print("  x, z: Spatial coordinates")

# Execute the function to print the result
print_magnetic_field_expression()