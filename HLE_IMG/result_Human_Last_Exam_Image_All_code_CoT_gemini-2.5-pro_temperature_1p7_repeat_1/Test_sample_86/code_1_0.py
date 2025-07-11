import math

def print_energy_spectrum():
    """
    Prints the derived formula for the energy spectrum of the perturbed harmonic oscillator
    based on the self-energy diagram.
    """
    
    # The final formula for the energy spectrum E_n includes several constants.
    # The numbers in the final equation are 1, 2, and 4.
    
    line1 = "E_n = (n + 1 / 2) * hbar * sqrt(omega_0**2 + (u * hbar) / (4 * m**2 * omega_0))"

    print("The energy spectrum E_n predicted by the self-energy diagram is given by the following equation:")
    print("This result is derived by calculating the corrected oscillator frequency based on the tadpole diagram correction.")
    print("\nThe final formula, showing each number, is:")
    print(line1)
    
    # Example of printing each number explicitly
    print("\nBreaking down the constants in the equation:")
    print("The energy levels are proportional to (n + 1/2), where 1 and 2 are present.")
    print("The frequency correction term has a denominator of 4 * m**2 * omega_0, where 4 and 2 are present.")


print_energy_spectrum()
