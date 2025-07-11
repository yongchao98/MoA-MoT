import math

def print_absorption_cross_section_equations():
    """
    This function prints the equations for the absorption cross-section σ(ω)
    for a chain of molecules under two different assumptions, based on
    first-order time-dependent perturbation theory for a Gaussian pulse.
    """

    # --- Introduction and Definitions ---
    print("This script formulates the equations for the absorption cross-section σ(ω).")
    print("-" * 70)
    print("Definitions of symbols used:")
    print("  σ(ω)        : Absorption cross-section as a function of laser frequency ω.")
    print("  ω           : Carrier frequency of the Gaussian laser pulse.")
    print("  τ           : Duration (width) of the Gaussian laser pulse.")
    print("  ħ           : Reduced Planck's constant.")
    print("  C_const     : Proportionality constant, C_const = 2*sqrt(π)*τ / (ħ*c*ε_0).")
    print("  exp(x)      : The exponential function e^x.")
    print("  sum_{i,f}   : A sum over initial occupied orbitals (i) and final unoccupied orbitals (f).")
    print("  E_i, E_f    : Energies of the initial and final molecular orbitals.")
    print("  |μ_fi|²     : Squared magnitude of the transition dipole moment between orbitals i and f.")
    print("  N           : Number of molecules in the chain.")
    print("  J           : Near-neighbor interaction energy (coupling constant).")
    print("-" * 70)

    # --- Case a) No interaction ---
    print("a) The interaction between molecules can be neglected.\n")
    print("The total cross-section is the sum of the cross-sections of N independent molecules.")
    print("Each molecule has the same set of possible transitions from occupied orbitals (i) to")
    print("unoccupied orbitals (f). The final equation is:\n")

    # Equation for Case (a)
    equation_a = (
        "σ_a(ω) = N * C_const * sum_{i,f} [ "
        "( (E_f - E_i) / ħ ) * |μ_fi|² * "
        "exp( -( ( (E_f - E_i) / ħ - ω )² * τ² ) )"
        " ]"
    )
    print(equation_a)
    print("\n" + "-" * 70)

    # --- Case b) Near-neighbor interaction ---
    print("b) The interaction between near-neighbors should be considered.\n")
    print("Interaction delocalizes the excitation into exciton states. For a simple model of a long")
    print("chain where each molecule has one relevant transition (i -> f), a selection rule emerges:")
    print("only the k=0 exciton state is optically active. Its transition energy is shifted by the")
    print("interaction, and its transition dipole moment is coherently enhanced. The final equation is:\n")

    # Equation for Case (b)
    # The transition energy is E_f - E_i + 2J.
    # The squared transition dipole moment is N * |μ_fi|².
    equation_b = (
        "σ_b(ω) = C_const * [ "
        "( (E_f - E_i + 2J) / ħ ) * ( N * |μ_fi|² ) * "
        "exp( -( ( (E_f - E_i + 2J) / ħ - ω )² * τ² ) )"
        " ]"
    )
    print(equation_b)
    print("-" * 70)


if __name__ == '__main__':
    print_absorption_cross_section_equations()
