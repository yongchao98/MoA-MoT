import numpy as np

def illustrate_partner_hamiltonian_spectra():
    """
    Explains and demonstrates the difference in spectra for a pair of
    supersymmetric partner Hamiltonians using the quantum harmonic oscillator example.
    """
    print("This script illustrates the spectral properties of two supersymmetric partner Hamiltonians.")
    print("Theory states their spectra can differ by at most one energy level.")
    print("We will use the quantum harmonic oscillator as a concrete example.")
    print("-" * 50)
    
    # We choose a superpotential related to the harmonic oscillator.
    # In the convention H_0 = L^+L and H_1 = LL^+ (with L=d/dx+W), the choice
    # W(x) = omega * x leads to the partner potentials V_0 = omega^2*x^2 - omega
    # and V_1 = omega^2*x^2 + omega. The constant alpha is absorbed into the potentials here.
    # The Hamiltonians are H_0 = -d^2/dx^2 + V_0 and H_1 = -d^2/dx^2 + V_1.
    
    # Let's set a value for the parameter omega.
    omega = 2.0
    num_levels = 6

    print(f"Example system: Partner Hamiltonians from the Harmonic Oscillator with omega = {omega}\n")
    
    # The spectrum of the standard Harmonic Oscillator H = -d^2/dx^2 + omega^2*x^2 is (2n+1)*omega.
    # The spectrum of H_0 = H - omega is E_n = (2n+1)*omega - omega = 2*n*omega
    spec_H0 = [2 * n * omega for n in range(num_levels)]
    
    # The spectrum of H_1 = H + omega is E_n = (2n+1)*omega + omega = (2n+2)*omega
    spec_H1 = [(2 * n + 2) * omega for n in range(num_levels)]

    print(f"The first {num_levels} theoretical energy levels for Hamiltonian H_0 are:")
    print("Formula: E_n = 2 * n * omega")
    for n in range(num_levels):
        print(f"  n={n}: 2 * {n} * {omega:.1f} = {spec_H0[n]:.1f}")
        
    print(f"\nThe first {num_levels} theoretical energy levels for Hamiltonian H_1 are:")
    print("Formula: E_n = (2 * n + 2) * omega")
    for n in range(num_levels):
        print(f"  n={n}: (2 * {n} + 2) * {omega:.1f} = {spec_H1[n]:.1f}")

    print("\n" + "-" * 50)
    print("Comparing the infinite sets of eigenvalues:")
    print("Spectrum of H_0: {0.0, 4.0, 8.0, 12.0, ...}")
    print("Spectrum of H_1: {4.0, 8.0, 12.0, 16.0, ...}")
    
    print("\nThe spectrum of H_0 contains the ground state energy level 0.0, which is absent from the spectrum of H_1.")
    print("All other energy levels of H_0 are present in H_1, and vice-versa.")
    print("Thus, the spectra differ by exactly one level.")
    print("\nThis confirms that the maximum number of differing levels is at least 1.")
    print("The theory guarantees that this number can't exceed 1.")


# Execute the function to show the results
illustrate_partner_hamiltonian_spectra()

print("\n<<<1>>>")