def solve_energy_spectrum():
    """
    This function prints the derived energy spectrum for the perturbed harmonic oscillator.
    """
    # Define the symbols as strings for clear printing
    E_n = "E_n"
    E_0 = "E_0"
    n = "n"
    hbar = "ħ"
    omega_0 = "ω₀"
    u = "u"
    m = "m"

    # Print the final equation for the energy spectrum
    print("The energy spectrum predicted by the self-energy diagram, up to an overall constant E₀, is:")
    # The final equation is E_n = E_0 + n*hbar*omega_prime
    # We substitute the calculated value for omega_prime
    print(f"{E_n} = {E_0} + {n} * {hbar} * {omega_0} + {n} * ({u} * {hbar}^3) / (8 * {m}^2 * {omega_0}^2)")

solve_energy_spectrum()