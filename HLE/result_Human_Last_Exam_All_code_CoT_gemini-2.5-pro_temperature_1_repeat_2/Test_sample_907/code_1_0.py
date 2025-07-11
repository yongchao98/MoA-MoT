def print_absorption_cross_section_equations():
    """
    This script prints the equations for the absorption cross-section of a molecular chain
    interacting with an ultrashort Gaussian laser pulse, based on first-order time-dependent
    perturbation theory.
    """

    print("----------------------------------------------------------------------")
    print("Equations for Absorption Cross-Section")
    print("----------------------------------------------------------------------\n")

    print("Symbols used in the equations:")
    print("  σ(ω) : Absorption cross-section as a function of laser frequency ω.")
    print("  C    : A constant containing fundamental constants (ħ, c, ε₀) and other factors.")
    print("  ω    : Angular frequency of the incident laser pulse.")
    print("  ω₀   : Transition angular frequency of an isolated molecule (HOMO-LUMO gap).")
    print("  τ    : Duration of the Gaussian laser pulse.")
    print("  N    : Number of molecules in the chain.")
    print("  μ    : Transition dipole moment of a single molecule.")
    print("  ε    : Polarization vector of the laser's electric field.")
    print("  J    : Nearest-neighbor coupling energy (interaction strength).")
    print("  ħ    : Reduced Planck's constant.")
    print("  j    : Quantum number for the exciton state (j = 1, 2, 3, ...).")
    print("\n----------------------------------------------------------------------\n")


    # --- Case a) No interaction between molecules ---
    print("a) Case where interaction between molecules is neglected:\n")
    print("The absorption spectrum consists of a single Gaussian peak. The cross-section is the")
    print("incoherent sum of the absorption from N identical molecules.\n")

    # Equation for case a
    equation_a = (
        "σ_a(ω) = C ⋅ ω ⋅ N ⋅ |μ ⋅ ε|² ⋅ exp( - (ω₀ - ω)² ⋅ τ² / 2 )"
    )
    print("Equation:\n")
    print(f"  {equation_a}\n")
    print("----------------------------------------------------------------------\n")


    # --- Case b) Nearest-neighbor interaction ---
    print("b) Case where nearest-neighbor interaction is considered:\n")
    print("The interaction leads to the formation of delocalized exciton states. The absorption")
    print("spectrum is a sum over the optically allowed exciton states (those with odd quantum")
    print("numbers 'j'), leading to multiple peaks within an 'exciton band'.\n")

    # Equation for case b
    equation_b_main = (
        "σ_b(ω) = C ⋅ ω ⋅ Σ_{j=1,3,...}^N |μ ⋅ ε|² ⋅ (2 / (N+1)) ⋅ cot²(πj / (2(N+1))) ⋅ exp( - (ω_j - ω)² ⋅ τ² / 2 )"
    )
    equation_b_energy = (
        "ω_j = ω₀ + (2J / ħ) ⋅ cos(πj / (N+1))"
    )
    print("Equation:\n")
    print(f"  {equation_b_main}\n")
    print("where the frequency of each exciton peak, ω_j, is given by:\n")
    print(f"  {equation_b_energy}\n")
    print("----------------------------------------------------------------------")


if __name__ == '__main__':
    print_absorption_cross_section_equations()