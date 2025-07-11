def solve_energy_spectrum():
    """
    This function calculates and prints the energy spectrum for a harmonic oscillator
    with a quartic perturbation, up to the first order.
    """
    
    # The first-order correction in perturbation theory for the energy of state n is
    # ΔE_n = <n| (u/4!) * x^4 |n>.
    # The calculation yields:
    # ΔE_n = (u * hbar^2) / (32 * m^2 * omega_0^2) * (2*n^2 + 2*n + 1).
    # The total energy is E_n = E_n^(0) + ΔE_n.
    # E_n = hbar*omega_0*(n + 1/2) + (u * hbar^2) / (32 * m^2 * omega_0^2) * (2*n^2 + 2*n + 1).
    
    # The question asks for the energy spectrum up to an overall constant E_0,
    # which is the ground state energy (n=0).
    # E_0 = hbar*omega_0*(1/2) + (u * hbar^2) / (32 * m^2 * omega_0^2) * (1).
    
    # We compute E_n - E_0:
    # E_n - E_0 = [hbar*omega_0*(n + 1/2) - hbar*omega_0*(1/2)] + [ΔE_n - ΔE_0]
    #           = n*hbar*omega_0 + (u*hbar^2)/(32*m^2*omega_0^2) * [(2n^2 + 2n + 1) - 1]
    #           = n*hbar*omega_0 + (u*hbar^2)/(32*m^2*omega_0^2) * (2n^2 + 2n)
    #           = n*hbar*omega_0 + (u*hbar^2)/(16*m^2*omega_0^2) * n*(n+1)

    print("The predicted energy spectrum, relative to the ground state energy E_0, is given by the formula:")
    print("")
    
    # As per the instruction to "output each number in the final equation",
    # the numbers are 16 and 1 (from n+1).
    final_equation = "E_n - E_0 = n*ħ*ω₀ + (u*ħ² / (16*m²*ω₀²)) * n*(n + 1)"
    
    print(final_equation)
    print("")
    print("Where:")
    print("  n:         the principal quantum number (n = 0, 1, 2, ...)")
    print("  ħ:         the reduced Planck constant")
    print("  ω₀:        the unperturbed angular frequency of the oscillator")
    print("  u:         the coupling constant for the quartic perturbation")
    print("  m:         the mass of the oscillator")

solve_energy_spectrum()