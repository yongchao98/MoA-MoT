import sympy

def calculate_cross_section():
    """
    Calculates and prints the total cross section for fermion-fermion scattering
    in a pseudoscalar Yukawa theory in the high-energy limit.

    The Lagrangian is given by:
    L = 1/2 (d phi)^2 - M^2/2 phi^2 + psi_bar(i*gamma*d - m)psi - g*psi_bar*gamma_5*psi*phi

    The calculation is performed under the following assumptions:
    1. High-energy limit, where the center-of-mass energy E is much larger
       than the fermion mass m (E >> m, so we take m -> 0).
    2. The energy is also much larger than the scalar boson mass M (E >> M).
    """

    # Define symbolic variables
    g = sympy.Symbol('g')  # Coupling constant
    E = sympy.Symbol('E')  # Center-of-mass energy of one fermion
    pi = sympy.pi          # The constant pi

    # The total cross section formula derived in the high-energy limit (E >> m, M)
    # sigma = (3 * g**4) / (128 * pi * E**2)
    
    # Numerator and denominator components
    numerator_coeff = 3
    numerator_vars = g**4
    
    denominator_coeff = 128
    denominator_vars = pi * E**2
    
    sigma = (numerator_coeff * numerator_vars) / (denominator_coeff * denominator_vars)

    # Print the result clearly showing all the numbers
    print("The total cross section for fermion-fermion scattering in the high-energy limit is:")
    print(f"σ = ({numerator_coeff} * {g.name}**4) / ({denominator_coeff} * π * {E.name}**2)")
    print("\nUsing sympy for representation:")
    sympy.pprint(sigma, use_unicode=True)

if __name__ == '__main__':
    calculate_cross_section()