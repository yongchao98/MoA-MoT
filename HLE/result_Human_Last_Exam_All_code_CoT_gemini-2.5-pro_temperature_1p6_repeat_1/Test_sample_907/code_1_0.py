import math

def print_absorption_equations():
    """
    This function prints the equations for the absorption cross-section
    for a chain of molecules under two different assumptions.
    """

    # --- Introduction and Definition of Symbols ---
    print("This script provides the equations for the absorption cross-section \u03C3(\u03C9) of a molecular chain.")
    print("The system absorbs an ultrashort Gaussian laser pulse of duration \u03C4 and carrier frequency \u03C9.")
    print("The derivation is based on first-order time-dependent perturbation theory.")
    print("\nThe symbols used in the equations are:")
    print("  \u03C3(\u03C9): Absorption cross-section as a function of laser frequency \u03C9.")
    print("  N: The number of molecules in the chain.")
    print("  \u03C9_eg: The transition frequency of a single, isolated molecule.")
    print("  \u03BC_eg: The transition dipole moment of a single, isolated molecule.")
    print("  \u03C4: The duration of the Gaussian laser pulse.")
    print("  J: The nearest-neighbor interaction energy (coupling constant).")
    print("  \u0127: The reduced Planck's constant.")
    print("  C: Proportionality constant grouping physical constants (\u221A\u03C0 / (2\u0127c\u03B5_0)).")
    print("-" * 60)

    # --- Case a) No interaction between molecules ---
    print("\nCase a): The interaction between molecules can be neglected.")
    print("The total absorption is the incoherent sum of absorptions by N identical molecules.")
    print("The peak absorption occurs at the single-molecule transition frequency \u03C9_eg.")
    print("\nThe equation for the absorption cross-section is:")

    # Equation format for case a)
    # The 'f' string allows embedding expressions inside a string.
    # Unicode characters are used for Greek letters.
    # Note: we use omega_eg inside the equation, but in the prefactor the frequency dependence
    # is often simplified to just omega, since it's a slow-varying function compared to the exponential.
    equation_a = "  \u03C3_a(\u03C9) = C \u22C5 N \u22C5 \u03C9 \u22C5 |\u03BC_eg|\u00B2 \u22C5 exp[ -(\u03C9 - \u03C9_eg)\u00B2 \u03C4\u00B2 ]"
    print(equation_a)
    print("-" * 60)

    # --- Case b) Nearest-neighbor interaction considered ---
    print("\nCase b): The interaction between near-neighbors should be considered.")
    print("The molecular excitations couple to form delocalized exciton states.")
    print("Optical selection rules allow transitions only to the k=0 exciton state.")
    print("This leads to a single absorption peak, shifted in energy by 2J, and a coherently enhanced transition strength.")
    print("\nThe equation for the absorption cross-section is:")

    # Equation format for case b)
    equation_b = "  \u03C3_b(\u03C9) = C \u22C5 N \u22C5 \u03C9 \u22C5 |\u03BC_eg|\u00B2 \u22C5 exp[ -(\u03C9 - (\u03C9_eg + 2J/\u0127))\u00B2 \u03C4\u00B2 ]"
    print(equation_b)


if __name__ == '__main__':
    print_absorption_equations()
