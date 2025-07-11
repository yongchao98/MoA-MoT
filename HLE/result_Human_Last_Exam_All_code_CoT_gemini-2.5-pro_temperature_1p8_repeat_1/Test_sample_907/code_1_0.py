import sys

def print_absorption_equations():
    """
    This function prints the derived equations for the absorption cross-section
    for a chain of molecules under two different conditions, as per the user's request.
    """

    # --- Introduction ---
    intro_text = """
The equation for the absorption cross-section, sigma(omega), is derived from
first-order time-dependent perturbation theory. For an excitation by an
ultrashort Gaussian laser pulse of duration tau, the otherwise sharp spectral
transitions are broadened. The resulting lineshape is a Gaussian function with
a spectral width inversely proportional to tau.

Below are the equations for the two requested cases.
"""
    print(intro_text)

    # --- Case a) No Interaction ---
    case_a_header = "--- Case a) Interaction between molecules can be neglected ---"
    case_a_explanation = """
In this case, all N molecules in the chain are independent. The total absorption
is the sum of the absorption of each molecule. All molecules have the same
transition frequency omega_eg. The total squared transition dipole moment is the
sum of the individual ones.
"""
    case_a_equation = "sigma_a(omega) = C * N * |mu_eg|^2 * omega * exp( -(tau**2 * (omega - omega_eg)**2) / 2 )"

    print(case_a_header)
    print(case_a_explanation)
    print("The equation for the absorption cross-section is:\n")
    print(f"  {case_a_equation}\n")

    # --- Case b) Nearest-Neighbor Interaction ---
    case_b_header = "--- Case b) Nearest-neighbor interaction should be considered ---"
    case_b_explanation = """
With nearest-neighbor interaction (coupling) J, the excitations become delocalized
Frenkel excitons. The new eigenstates form an energy band. Due to selection rules,
only the exciton state with wavevector k=0 is optically active (bright). This state's
energy is shifted by the coupling J. The transition strength to this single state
is coherently enhanced by a factor of N.
"""
    # Using 'hbar' for the reduced Planck constant.
    case_b_equation = "sigma_b(omega) = C * N * |mu_eg|^2 * omega * exp( -(tau**2 * (omega - (omega_eg + 2*J/hbar))**2) / 2 )"

    print(case_b_header)
    print(case_b_explanation)
    print("The equation for the absorption cross-section is:\n")
    print(f"  {case_b_equation}\n")

    # --- Symbol Definitions ---
    definitions = """
Where the symbols in the equations represent:
  sigma(omega) : Absorption cross-section as a function of light frequency.
  C            : A proportionality constant that includes fundamental physical constants.
  N            : The number of molecules in the chain.
  |mu_eg|^2    : The squared transition dipole moment for a single molecule's
                 electronic transition (from ground 'g' to excited 'e' state).
  omega        : The angular frequency of the absorbed light.
  tau          : The duration of the ultrashort Gaussian laser pulse.
  omega_eg     : The transition angular frequency for an isolated molecule.
  J            : The nearest-neighbor coupling energy.
  hbar         : The reduced Planck constant.
  exp()        : The exponential function, which defines the Gaussian lineshape.
"""
    print(definitions)
    
    # Final answer tag as per instructions
    sys.stdout.write("<<<")
    sys.stdout.write(f"Equation for a): {case_a_equation}. Equation for b): {case_b_equation}.")
    sys.stdout.write(">>>")


if __name__ == '__main__':
    print_absorption_equations()
