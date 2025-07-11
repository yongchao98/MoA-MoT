import math

def print_energy_spectrum():
    """
    This function prints the energy spectrum for a harmonic oscillator with a quartic perturbation,
    as predicted by the first-order self-energy tadpole diagram.
    """
    
    # The result is symbolic. We will represent it using strings.
    # The key result is the renormalized frequency, let's call it omega_prime.
    
    # The square of the renormalized frequency
    # (omega_prime)^2 = (omega_0)^2 + (u * hbar) / (4 * m^2 * omega_0)
    
    print("Based on the self-energy diagram, the system is equivalent to a new harmonic oscillator with a renormalized frequency, omega'.")
    print("-" * 80)
    
    # Define the numerator and denominator of the correction term for clarity
    numerator = "u * hbar"
    # Note the number 4 is explicitly included as requested.
    denominator = "({} * m**2 * omega_0)".format(4)
    
    print("The renormalized frequency-squared (omega')^2 is given by:")
    print("    (omega')^2 = (omega_0)^2 + ({}) / {}".format(numerator, denominator))
    print()

    # The energy spectrum is E_n = E_0 + n * hbar * omega'
    print("The energy spectrum E_n, up to an overall constant E_0, is predicted to be:")
    print("    E_n = E_0 + n * hbar * sqrt((omega_0)^2 + ({}) / {})".format(numerator, denominator))
    print("-" * 80)
    print("Where:")
    print("    E_n: Energy of the n-th level")
    print("    E_0: Constant ground state energy")
    print("    n: Quantum number (0, 1, 2, ...)")
    print("    hbar: Reduced Planck's constant")
    print("    omega_0: Original angular frequency")
    print("    u: Perturbation strength")
    print("    m: Mass of the oscillator")

if __name__ == "__main__":
    print_energy_spectrum()
