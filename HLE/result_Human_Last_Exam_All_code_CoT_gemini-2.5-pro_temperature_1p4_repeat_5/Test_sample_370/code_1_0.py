import sympy

def calculate_cross_section():
    """
    Calculates and prints the total cross-section for fermion-fermion scattering
    in the specified theory, in the high-energy limit.
    """
    # Define symbols for the variables in the equation
    g, E, pi = sympy.symbols('g E pi')

    # The formula for the total cross section sigma is derived as:
    # sigma = (3 * g^4) / (128 * pi * E^2)
    
    # Define the numerical coefficients from the derived formula
    numerator_coeff = 3
    denominator_coeff = 128
    
    # Construct the symbolic expression for the cross-section
    sigma = (numerator_coeff * g**4) / (denominator_coeff * pi * E**2)
    
    print("The total cross section (sigma) for the scattering of two fermions in the high-energy limit is:")
    
    # Use sympy.pretty_print for a nicely formatted output of the equation
    sympy.pprint(sympy.Eq(sympy.Symbol('sigma'), sigma), use_unicode=True)
    
    print("\nIn this final equation:")
    print(f"The number in the numerator is: {numerator_coeff}")
    print(f"The number in the denominator is: {denominator_coeff}")
    print("It depends on the coupling constant 'g' to the power of 4.")
    print("It is inversely proportional to the square of the energy 'E'.")

if __name__ == '__main__':
    calculate_cross_section()
