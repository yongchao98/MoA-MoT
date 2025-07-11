import math

def print_force_formula():
    """
    This function prints the components of the derived formula for the force per unit area.
    The formula is symbolic, representing the physical quantities involved.
    """
    
    print("Based on the derivation, the force per unit area, vec(f), on the plane at x=d is calculated.")
    print("Direction: The force is in the positive x-direction (i_x).")
    print("-" * 30)
    print("Magnitude Equation: f = Numerator / Denominator")
    print("-" * 30)
    
    # Describing the numerator
    print("Numerator = (1/2) * mu_0 * K_0^2 * cos^2(omega * t)")
    print("  - Factor: 1/2")
    print("    - Number '1'")
    print("    - Number '2'")
    print("  - mu_0: Permeability of free space")
    print("  - K_0^2: Surface current amplitude squared")
    print("  - cos^2(omega * t): Time-varying component with power '2'")
    
    print("-" * 30)
    
    # Describing the denominator
    print("Denominator = cosh^2(omega_p * d / c)")
    print("  - cosh^2: Hyperbolic cosine squared, with power '2'")
    print("  - Argument: (omega_p * d / c)")
    
    print("-" * 30)
    print("Final Derived Formula:")
    print("vec(f) = i_x * (1/2) * mu_0 * K_0^2 * cos^2(omega * t) / cosh^2(omega_p * d / c)")

# Execute the function to display the formula components
print_force_formula()