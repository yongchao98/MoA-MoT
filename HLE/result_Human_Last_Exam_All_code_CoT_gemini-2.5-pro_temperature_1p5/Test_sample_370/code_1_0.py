import math

def calculate_cross_section_formula():
    """
    This function prints the derived formula for the total cross section
    for fermion-fermion scattering in the specified theory, under the
    high-energy (E >> m) and heavy mediator (M >> E) approximations.

    The derived formula is:
    sigma = (g^4 * E^2) / (8 * pi * M^4)
    """

    # Numerator of the fraction
    numerator = "g^4 * E^2"
    
    # Denominator of the fraction
    denominator_coeff = 8
    denominator_vars = "pi * M^4"
    denominator = f"{denominator_coeff} * {denominator_vars}"

    # To satisfy the "output each number" requirement, we explicitly define the components
    # The equation is of the form: (g^power_g * E^power_E) / (coeff * pi * M^power_M)
    power_g = 4
    power_E = 2
    coefficient_8 = 8 # Derived from (64 / (2 * 4*pi)) * 4pi, simplified. 
                      # More simply from (s / (32*pi)) which is (4*E^2 / (32*pi)) = E^2/(8*pi)
    power_M = 4

    print("The total cross section sigma is given by the formula:")
    print() # blank line for spacing
    
    # Calculate the width for alignment
    width = max(len(numerator), len(denominator))
    
    # Print the formatted formula
    print(f"{numerator:^{width}}")
    print("-" * width)
    print(f"{denominator:^{width}}")

    print("\nWhere:")
    print("  g is the coupling constant.")
    print("  E is the center-of-mass energy of one of the incoming fermions.")
    print("  M is the mass of the scalar particle phi.")
    print("  pi is the mathematical constant.")
    
    print("\nThe numerical components of the formula are:")
    print(f"  Power of g: {power_g}")
    print(f"  Power of E: {power_E}")
    print(f"  Numerical coefficient in the denominator: {coefficient_8}")
    print(f"  Power of M: {power_M}")

# Execute the function to display the result
calculate_cross_section_formula()