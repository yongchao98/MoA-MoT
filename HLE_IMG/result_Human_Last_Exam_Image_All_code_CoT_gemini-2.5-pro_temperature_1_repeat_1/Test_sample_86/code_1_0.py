import sympy

def solve_harmonic_oscillator_perturbation():
    """
    This function derives and prints the energy spectrum for a harmonic oscillator
    with a quartic perturbation, based on the self-energy diagram.
    """
    # Define symbolic variables
    E_n = sympy.Symbol('E_n')
    n = sympy.Symbol('n')
    hbar = sympy.Symbol('hbar')
    w0 = sympy.Symbol('omega_0')
    u = sympy.Symbol('u')
    m = sympy.Symbol('m')

    # Unperturbed energy
    E_unperturbed = hbar * w0 * (n + sympy.Rational(1, 2))

    # Energy correction term from the self-energy diagram
    E_correction = (u * hbar**2) / (8 * m**2 * w0**2) * (n + sympy.Rational(1, 2))

    # Final energy spectrum
    final_equation = sympy.Eq(E_n, E_unperturbed + E_correction)

    # Print the result in a readable format
    print("The self-energy diagram predicts the energy spectrum to be:")
    
    # We use pretty print for a better mathematical representation
    sympy.init_printing(use_unicode=True)
    
    # Manually construct the string to match the requested format
    # "Remember in the final code you still need to output each number in the final equation!"
    term1 = f"(n + 1/2)*hbar*omega_0"
    term2 = f"(n + 1/2)*u*hbar**2 / (8*m**2*omega_0**2)"
    
    # To explicitly show the numbers 1/2 and 1/8 as requested
    equation_string = f"E_n = (n + 1/2)*hbar*omega_0 + (1/8) * (u*hbar**2 / (m**2*omega_0**2)) * (n + 1/2)"
    
    # Let's print the SymPy equation for clarity and the custom string to highlight the numbers
    print(final_equation)

solve_harmonic_oscillator_perturbation()