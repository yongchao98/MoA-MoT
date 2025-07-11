import sympy as sp

def display_absorption_equations():
    """
    Generates and prints the equations for the absorption cross-section
    for a chain of molecules absorbing a Gaussian laser pulse, under
    two different theoretical assumptions.
    """

    # --- Define all the physical symbols used in the equations ---
    sigma = sp.Symbol('σ(ω)', real=True, positive=True)
    sigma_agg = sp.Symbol('σ_agg(ω)', real=True, positive=True)
    omega = sp.Symbol('ω', real=True, positive=True)
    omega_fi = sp.Symbol('ω_fi', real=True, positive=True)
    tau = sp.Symbol('τ', real=True, positive=True)
    d_fi_sq = sp.Symbol('|d_fi|²', real=True, positive=True)
    epsilon_0 = sp.Symbol('ε₀', real=True, positive=True)
    c = sp.Symbol('c', real=True, positive=True)
    hbar = sp.Symbol('ħ', real=True, positive=True)
    J = sp.Symbol('J', real=True)
    N = sp.Symbol('N', real=True, integer=True, positive=True)

    # Use SymPy's pretty printing for clear mathematical output
    sp.init_printing(use_unicode=True)

    # --- Part a) No Interaction Between Molecules ---
    print("="*60)
    print("Part a) The interaction between molecules can be neglected.")
    print("="*60)
    print("In this scenario, we consider each molecule independently. The absorption")
    print("cross-section, σ(ω), describes the probability of a transition from the")
    print("ground state to an excited state in a single molecule, induced by a")
    print("Gaussian laser pulse. The shape of the absorption spectrum is determined")
    print("by the Fourier transform of the pulse.")
    print("\nThe equation for the absorption cross-section is:")

    # Construct the equation for a single molecule
    prefactor_a = (sp.sqrt(sp.pi) * omega * tau * d_fi_sq) / (epsilon_0 * c * hbar)
    exponent_a = -(omega - omega_fi)**2 * tau**2
    equation_a = sp.Eq(sigma, prefactor_a * sp.exp(exponent_a))

    print()
    print(sp.pretty(equation_a))
    print()

    print("Where each symbol in the equation represents:")
    print(f"  σ(ω)      : Absorption cross-section as a function of light frequency ω")
    print(f"  ω         : Frequency of the incident laser light")
    print(f"  ω_fi      : Transition frequency of a single molecule, (E_final - E_initial)/ħ")
    print(f"  τ         : Temporal duration (width) of the Gaussian laser pulse")
    print(f"  |d_fi|²   : Squared magnitude of the transition dipole moment of one molecule")
    print(f"  ε₀        : Permittivity of free space")
    print(f"  c         : Speed of light in vacuum")
    print(f"  ħ         : Reduced Planck constant")

    # --- Part b) Nearest-Neighbor Interaction Considered ---
    print("\n\n" + "="*60)
    print("Part b) The interaction between near-neighbors should be considered.")
    print("="*60)
    print("When nearest-neighbor interactions are included, the excitations are")
    print("delocalized over the molecular chain, forming 'Frenkel excitons'. A key")
    print("selection rule arises: only the collective exciton state with wavevector k=0")
    print("is optically bright. This changes the absorption properties significantly.")
    print("\nThe equation for the aggregate absorption cross-section is:")

    # Construct the equation for the interacting chain (aggregate)
    # The transition frequency is shifted by the coupling energy
    shifted_omega = omega_fi + (2 * J / hbar)
    # The oscillator strength is enhanced by the number of molecules
    prefactor_b = (sp.sqrt(sp.pi) * omega * tau * N * d_fi_sq) / (epsilon_0 * c * hbar)
    exponent_b = -(omega - shifted_omega)**2 * tau**2
    equation_b = sp.Eq(sigma_agg, prefactor_b * sp.exp(exponent_b))

    print()
    print(sp.pretty(equation_b))
    print()

    print("Where each new symbol in the equation represents:")
    print(f"  σ_agg(ω)  : Absorption cross-section of the entire molecular aggregate")
    print(f"  N         : The number of molecules coherently coupled in the chain")
    print(f"  J         : The nearest-neighbor coupling energy")
    print("\nNote the key differences from the non-interacting case:")
    print(f"1. Absorption strength is enhanced by a factor of N (superradiance).")
    print(f"2. The absorption peak is shifted in frequency by 2J/ħ.")

if __name__ == "__main__":
    display_absorption_equations()