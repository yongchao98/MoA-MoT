def print_absorption_cross_section_formulas():
    """
    This script prints the equations for the absorption cross-section of a molecular chain
    interacting with an ultrashort Gaussian laser pulse, derived from first-order
    time-dependent perturbation theory.
    """

    print("="*80)
    print("Absorption Cross-Section for a Molecular Chain")
    print("="*80)
    print("\nThe general form for the absorption cross-section sigma(omega) is derived from")
    print("first-order perturbation theory for a Gaussian pulse. It is given by:")
    print("\nsigma(omega) = C * omega * sum_{i,f} |mu_fi|^2 * exp(-(omega_fi - omega)^2 * tau^2)\n")
    print("We will now specify this general form for the two requested cases.")

    # Case a: No interaction between molecules
    print_case_a()

    # Case b: Near-neighbor interaction
    print_case_b()

def print_case_a():
    """
    Prints the formula and symbol definitions for case (a):
    Interaction between molecules can be neglected.
    """
    equation_a = (
        "sigma_a(omega) = (2 * pi^(3/2) * tau) / (3 * hbar * c * epsilon_0) * omega * "
        "Sum_{j,l} |mu_lj|^2 * exp(-((epsilon_l - epsilon_j) / hbar - omega)^2 * tau^2)"
    )

    print("\n" + "-"*80)
    print("a) Case: Interaction between molecules is neglected.")
    print("-"*80)
    print("\nIn this scenario, the total absorption is the sum of absorptions from individual,")
    print("identical molecules. The energy levels are the discrete molecular orbitals.")
    print("\nThe equation for the absorption cross-section is:")
    print(f"\n{equation_a}\n")

    print("Where the symbols represent:")
    print("  sigma_a(omega) : Absorption cross-section as a function of laser frequency omega.")
    print("  omega          : The central angular frequency of the laser pulse.")
    print("  pi             : The mathematical constant pi (approx. 3.14159).")
    print("  tau            : The duration of the Gaussian laser pulse.")
    print("  hbar           : The reduced Planck constant.")
    print("  c              : The speed of light in vacuum.")
    print("  epsilon_0      : The vacuum permittivity.")
    print("  Sum_{j,l}      : A sum over all occupied molecular orbitals (indexed by j)")
    print("                   and all unoccupied molecular orbitals (indexed by l).")
    print("  epsilon_j, epsilon_l : The energies of the respective molecular orbitals j and l.")
    print("  mu_lj          : The transition dipole moment between orbitals j and l.")
    print("  exp()          : The exponential function, which describes the Gaussian shape")
    print("                   of the absorption peaks due to the finite pulse duration.")


def print_case_b():
    """
    Prints the formula and symbol definitions for case (b):
    Near-neighbor interaction is considered.
    """
    equation_b = (
        "sigma_b(omega) = (2 * pi^(3/2) * tau) / (3 * hbar * c * epsilon_0) * omega * |mu_cv|^2 * "
        "Sum_k exp(-(((E_c - E_v - 2*(J_c - J_v)*cos(k*a)) / hbar) - omega)^2 * tau^2)"
    )
    print("\n" + "-"*80)
    print("b) Case: Interaction between near-neighbors is considered.")
    print("-"*80)
    print("\nHere, molecular orbitals hybridize to form delocalized bands (exciton bands).")
    print("Transitions occur between these bands, and the absorption spectrum becomes broad.")
    print("\nThe equation for a transition between a single valence and conduction band is:")
    print(f"\n{equation_b}\n")

    print("Where the symbols represent:")
    print("  sigma_b(omega) : Absorption cross-section as a function of laser frequency omega.")
    print("  omega, pi, tau, hbar, c, epsilon_0: Same as in case (a).")
    print("  mu_cv          : The transition dipole moment of a single molecule between the")
    print("                   occupied (valence, v) and unoccupied (conduction, c) orbitals")
    print("                   that form the bands.")
    print("  Sum_k          : A sum over all allowed wavevectors k in the first Brillouin zone")
    print("                   of the molecular chain.")
    print("  E_c, E_v       : The on-site energies of the molecular orbitals that form the")
    print("                   conduction and valence bands, respectively.")
    print("  J_c, J_v       : The near-neighbor coupling (transfer) integrals for the")
    print("                   conduction and valence bands.")
    print("  k              : The wavevector, which labels the states within a band.")
    print("  a              : The distance between adjacent molecules in the chain.")
    print("  cos()          : The cosine function, which describes the energy dispersion of the bands.")
    print("  exp()          : The exponential function, representing the Gaussian contribution")
    print("                   from the finite pulse duration.")


if __name__ == "__main__":
    print_absorption_cross_section_formulas()