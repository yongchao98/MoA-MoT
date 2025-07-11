import math

def calculate_energy_shift():
    """
    This script calculates and displays the ground state energy shift (van der Waals interaction)
    for two distant, identical quantum harmonic oscillators.

    The derivation uses second-order perturbation theory on the dipole-dipole interaction
    potential. The model assumes two 3D isotropic oscillators.
    """

    print("Derivation of the ground state energy shift (London Dispersion Force):")
    print("1. The interaction potential is the dipole-dipole interaction, H'.")
    print("2. The first-order energy correction <0|H'|0> is zero.")
    print("3. The second-order correction ΔE = Σ |<n|H'|0>|^2 / (E₀ - Eₙ) is calculated.")
    print("4. Summing over all contributing excited states for 3D oscillators gives the final result.")
    print("\nThe final formula for the leading term of the ground state energy shift ΔE is:")

    # The formula is: ΔE = - (C_num * e^4 * hbar) / (C_den * pi^2 * m^2 * omega_0^3 * R^6)
    # From the derivation, the factor is -6 * (1 / (2*hbar*omega_0)) * (1 / (16*pi^2)) * (hbar^2 / (4*m^2*omega_0^2))
    # = - (6 * hbar) / (128 * pi^2 * m^2 * omega_0^3)
    # = - (3 * hbar) / (64 * pi^2 * m^2 * omega_0^3)

    numerator_coeff = 3
    denominator_coeff = 64

    # Define variables with unicode characters for a clean mathematical expression
    e = "e"
    hbar = "ħ"
    pi = "π"
    m = "m"
    omega_0 = "ω₀"
    R = "R"

    # Construct and print the final equation string
    equation = f"ΔE = - ({numerator_coeff} * {e}⁴ * {hbar}) / ({denominator_coeff} * {pi}² * {m}² * {omega_0}³ * {R}⁶)"
    print(equation)

# Execute the function
calculate_energy_shift()