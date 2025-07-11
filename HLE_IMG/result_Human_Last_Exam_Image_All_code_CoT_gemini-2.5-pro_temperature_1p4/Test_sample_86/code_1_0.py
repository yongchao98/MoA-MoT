import sympy as sp

def solve_harmonic_oscillator_spectrum():
    """
    This function derives and prints the energy spectrum for a harmonic oscillator
    with a quartic perturbation, based on the self-energy diagram.
    """
    # Define symbols for the variables in the equation
    n = sp.Symbol('n')
    hbar = sp.Symbol('hbar')
    omega_0 = sp.Symbol('omega_0')
    u = sp.Symbol('u')
    m = sp.Symbol('m')
    E_n = sp.Symbol('E_n')

    # Define the numerical constants present in the formula
    num_1 = 1
    num_2 = 2
    num_4 = 4

    # Construct the string for the final equation as per the derivation
    # E_n = (n + 1/2) * hbar * sqrt(omega_0^2 + (u*hbar) / (4*m^2*omega_0))
    equation_string = (
        f"E_n = (n + {num_1}/{num_2}) * hbar * "
        f"sqrt(omega_0**2 + (u * hbar) / ({num_4} * m**2 * omega_0))"
    )

    print("The predicted energy spectrum up to an overall constant is:")
    print(equation_string)

solve_harmonic_oscillator_spectrum()