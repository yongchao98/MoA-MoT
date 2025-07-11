def solve_absorption_equation():
    """
    This script presents the equations for the absorption cross-section of a molecular chain
    interacting with an ultrashort Gaussian laser pulse, based on first-order time-dependent
    perturbation theory.
    """

    print("--- Absorption Cross-Section Equations ---")
    print("\nBelow are the equations for the two requested cases. The 'ultrashort Gaussian pulse' implies")
    print("that the absorption lineshape is Gaussian, with a width determined by the pulse duration τ.")
    print("-" * 40)

    # Case a) No interaction
    print("\na) Case a: No interaction between molecules.")
    print("   In this case, the chain is a simple collection of N independent molecules.")
    print("   The total absorption cross-section is N times that of a single molecule, and the spectrum")
    print("   shows a single absorption peak at the molecular transition frequency.")
    print("\n   The equation is:")
    print("   σ(ω) = N * (π * ω / (ε₀ * c)) * |μ_eg ⋅ ε|² * G(ω, ω_eg)")
    print("\n   Where:")
    print("     - σ(ω): Absorption cross-section as a function of light frequency ω.")
    print("     - N:     Number of molecules in the chain.")
    print("     - π:     The number Pi (3.14159...).")
    print("     - ε₀:    Permittivity of free space.")
    print("     - c:     Speed of light.")
    print("     - μ_eg:  Transition dipole moment for the electronic transition in a single molecule.")
    print("     - ε:     Polarization vector of the laser electric field.")
    print("     - ω_eg:  Angular frequency of the transition (E_eg / ħ).")
    print("     - G(ω, ω_eg): The normalized Gaussian lineshape from the pulse, defined as:")
    print("          G(ω, ω_eg) = (τ / sqrt(π)) * exp[-τ² * (ω - ω_eg)²]")
    print("     - τ:     A measure of the laser pulse duration.")

    print("-" * 40)

    # Case b) Near-neighbor interaction
    print("\nb) Case b: Near-neighbor interaction is considered (Frenkel Exciton Model).")
    print("   In this case, the interaction (coupling J) delocalizes the excitation into a band of")
    print("   N exciton states. Optical selection rules allow transitions only to specific exciton states.")
    print("\n   The equation is a sum over the allowed exciton transitions (k=1, 3, 5,...):")
    print("   σ(ω) = (π * ω / (ε₀ * c)) * |μ_eg ⋅ ε|² * Σ_k [F_k * G(ω, ω_k)]")
    print("\n   Where:")
    print("     - The sum Σ_k is over odd integer quantum numbers: k = 1, 3, 5, ... up to N.")
    print("     - F_k: The dimensionless oscillator strength for the transition to exciton state k.")
    print("          F_k = (2 / (N + 1)) * cot²(π * k / (2 * (N + 1)))")
    print("          Note: This term is largest for k=1, making it the dominant peak.")
    print("     - G(ω, ω_k): The Gaussian lineshape for the transition to exciton state k:")
    print("          G(ω, ω_k) = (τ / sqrt(π)) * exp[-τ² * (ω - ω_k)²]")
    print("     - ω_k: The angular frequency of the k-th exciton transition (E_k / ħ).")
    print("     - E_k: The energy of the k-th exciton state, given by:")
    print("          E_k = E_eg + 2 * J * cos(π * k / (N + 1))")
    print("     - J: The coupling energy between neighboring molecules.")
    print("     - E_eg: The transition energy of an isolated molecule.")

if __name__ == '__main__':
    solve_absorption_equation()