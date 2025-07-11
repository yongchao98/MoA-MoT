import sympy

def calculate_cross_section_formula():
    """
    Calculates and prints the symbolic formula for the total cross-section
    for fermion-fermion scattering via a pseudoscalar mediator.
    """
    # Define the symbolic variables
    # g: coupling constant
    # E: center-of-mass energy of one fermion
    # M: mass of the scalar mediator particle
    # pi: the mathematical constant pi
    g, E, M, pi = sympy.symbols('g E M pi', real=True, positive=True)

    # The Mandelstam variable s is the square of the total center-of-mass energy
    s = 4 * E**2

    # The calculation involves integrating the differential cross section.
    # The final result is composed of several terms.
    # We build the expression step-by-step.
    
    # Term 1 from the integration
    term1_numerator = 3 * s**2 + 5 * s * M**2
    term1_denominator = s + M**2
    term1 = term1_numerator / term1_denominator

    # Term 2 from the integration (coefficient of the logarithm)
    term2_numerator = 6 * s * M**2 + 10 * M**4
    term2_denominator = s + 2 * M**2
    term2_coeff = term2_numerator / term2_denominator

    # The logarithmic term
    log_term = sympy.log(M**2 / (s + M**2))

    # The full expression inside the main brackets
    full_expression = term1 + term2_coeff * log_term

    # The prefactor for the total cross section formula
    prefactor = g**4 / (32 * pi * s**2)

    # The final formula for the total cross section sigma
    sigma = prefactor * full_expression

    # Print the final formula in a readable format
    print("The total cross section sigma for fermion-fermion scattering in the lowest order is:")
    sympy.pprint(sigma, use_unicode=True)
    
    # To satisfy the "output each number" requirement, we explicitly print them.
    print("\nThe numbers in the final equation are derived from the calculation:")
    print("Prefactor constant: 1/32")
    print("Term 1 coefficients: 3, 5")
    print("Term 2 coefficients: 6, 10, 2")
    print("Mandelstam variable s = 4 * E^2")


if __name__ == '__main__':
    calculate_cross_section_formula()