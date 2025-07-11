def print_energy_spectrum():
    """
    Prints the derived energy spectrum of the perturbed harmonic oscillator.
    The formula is printed symbolically.
    """

    # The energy spectrum E_n consists of the original energy plus a correction term.
    # E_n = E_n^(0) + Delta_E_n
    
    # The original unperturbed energy levels are:
    # E_n^(0) = (n + 1/2) * hbar * omega_0
    
    # The correction due to the quartic perturbation, as predicted by the
    # self-energy diagram in the Hartree-Fock approximation is:
    # Delta_E_n = (n + 1/2) * (u * hbar^2) / (8 * m^2 * omega_0^2)

    print("The predicted energy spectrum E_n, up to an overall constant, is given by:")
    print("\nE_n = (n + 1/2) * hbar * omega_0 + (n + 1/2) * (u * hbar^2) / (8 * m^2 * omega_0^2)\n")

    print("Let's break down the final equation term by term as requested:")
    print("Term 1 (Unperturbed Energy): (n + 1/2) * hbar * omega_0")
    print(" - n: The energy level number (an integer >= 0)")
    print(" - 1/2: The zero-point energy factor, which is the number 0.5")
    print(" - hbar: The reduced Planck constant")
    print(" - omega_0: The natural frequency of the unperturbed oscillator")

    print("\nTerm 2 (Energy Correction): (n + 1/2) * (u * hbar^2) / (8 * m^2 * omega_0^2)")
    print(" - (n + 1/2): The same factor as above")
    print(" - u: The coupling constant of the quartic perturbation")
    print(" - hbar^2: The reduced Planck constant squared")
    print(" - 8: A numerical factor from the calculation, which is the number 8")
    print(" - m^2: The mass of the oscillator squared")
    print(" - omega_0^2: The unperturbed natural frequency squared")

if __name__ == '__main__':
    print_energy_spectrum()
