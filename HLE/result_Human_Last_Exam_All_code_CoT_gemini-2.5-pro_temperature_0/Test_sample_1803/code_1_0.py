import sympy

def calculate_energy_shift():
    """
    Calculates and prints the ground state energy shift for two interacting
    quantum harmonic oscillators using symbolic mathematics.
    """
    # Define the symbolic variables
    e = sympy.Symbol('e')          # Electron charge
    R = sympy.Symbol('R')          # Distance between oscillators
    m = sympy.Symbol('m')          # Mass of the oscillator
    omega0 = sympy.Symbol('omega0')  # Angular frequency of the oscillator
    hbar = sympy.Symbol('hbar')      # Reduced Planck constant
    pi = sympy.Symbol('pi')        # Pi

    # The final formula for the energy shift is derived from second-order
    # perturbation theory.
    # Delta_E = - (3 * e**4 * hbar) / (64 * pi**2 * R**6 * m**2 * omega0**3)

    # Construct the numerator and denominator of the expression
    numerator_coeff = -3
    denominator_coeff = 64

    numerator_symbols = e**4 * hbar
    denominator_symbols = pi**2 * R**6 * m**2 * omega0**3

    # Create the full expression
    Delta_E = (numerator_coeff * numerator_symbols) / (denominator_coeff * denominator_symbols)

    # Print the final equation, showing each numerical part
    print("The ground state energy shift (Delta_E) is given by the equation:")
    
    # Using sympy.pretty_print for a more readable output
    print("\nDelta_E =")
    sympy.pretty_print(Delta_E)

    # Explicitly state the numbers in the final equation as requested
    print("\nBreaking down the final equation:")
    print(f"Numerator coefficient: {numerator_coeff}")
    print(f"Power of charge 'e': 4")
    print(f"Power of 'hbar': 1")
    print(f"Denominator coefficient: {denominator_coeff}")
    print(f"Power of 'pi': 2")
    print(f"Power of distance 'R': 6")
    print(f"Power of mass 'm': 2")
    print(f"Power of frequency 'omega0': 3")


if __name__ == '__main__':
    calculate_energy_shift()
    # The final expression is symbolic.
    # To provide a single answer as requested by the format,
    # we present the symbolic formula itself.
    final_answer = "-3*e**4*hbar/(64*pi**2*R**6*m**2*omega0**3)"
    # print(f"\n<<<{final_answer}>>>")