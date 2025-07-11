import sympy

def solve_absorption_cross_section():
    """
    This function generates and prints the equations for absorption cross-section
    for a chain of molecules under two conditions: non-interacting and interacting.
    """

    # --- Define Symbolic Variables ---
    # These symbols represent the physical quantities involved in the equations.
    C = sympy.Symbol('C')                     # Represents all proportionality constants
    N = sympy.Symbol('N', positive=True)      # Number of molecules in the chain
    omega = sympy.Symbol('omega', real=True)  # Central frequency of the laser pulse
    omega_0 = sympy.Symbol('omega_0', positive=True) # Excitation frequency of a single molecule
    tau = sympy.Symbol('tau', positive=True)      # Duration of the Gaussian pulse
    mu_eg_sq = sympy.Symbol('|μ_eg|^2', positive=True) # Square of the transition dipole moment for one molecule
    J = sympy.Symbol('J', real=True)          # Near-neighbor coupling energy
    hbar = sympy.Symbol('ħ', positive=True)    # Reduced Planck constant

    # --- Case a) No interaction between molecules ---

    print("a) Interaction between molecules can be neglected.")
    print("=" * 50)
    print("The total absorption is the sum of absorptions from N independent molecules.")
    print("The absorption peak is at the single-molecule transition frequency ω₀.")
    print("The final equation for the absorption cross-section σ_a(ω) is:")

    # Construct the equation for case (a)
    # The cross-section is proportional to N, ω₀, |μ_eg|², and a Gaussian lineshape.
    sigma_a = C * N * omega_0 * mu_eg_sq * sympy.exp(-(omega - omega_0)**2 * tau**2)
    equation_a = sympy.Eq(sympy.Symbol('σ_a(ω)'), sigma_a)

    # Print the beautifully formatted equation
    sympy.pprint(equation_a, use_unicode=True)
    print("\nWhere the symbols represent:")
    print(f"  σ_a(ω)    : Absorption cross-section")
    print(f"  C         : Proportionality constant")
    print(f"  N         : Number of molecules")
    print(f"  ω         : Laser frequency")
    print(f"  ω₀        : Molecular transition frequency")
    print(f"  |μ_eg|²   : Squared transition dipole moment of a single molecule")
    print(f"  τ         : Gaussian pulse duration")
    print("\n" * 2)


    # --- Case b) Near-neighbor interaction is considered ---

    print("b) Interaction between near-neighbors should be considered.")
    print("=" * 50)
    print("Interactions create delocalized excitons. A selection rule allows transitions only to the k=0 exciton state.")
    print("This transition is coherently enhanced and its frequency is shifted by the coupling energy J.")

    # Define the new, shifted exciton frequency ω_J
    omega_J = sympy.Symbol('ω_J')
    equation_omega_J = sympy.Eq(omega_J, omega_0 + (2 * J / hbar))
    print("\nThe transition occurs at a new frequency, ω_J, defined as:")
    sympy.pprint(equation_omega_J, use_unicode=True)

    print("\nThe final equation for the absorption cross-section σ_b(ω) is:")
    # The oscillator strength |μ|² becomes |sqrt(N) * μ_eg|² = N * |μ_eg|²
    # The transition occurs at the shifted frequency ω_J.
    sigma_b = C * N * omega_J * mu_eg_sq * sympy.exp(-(omega - omega_J)**2 * tau**2)
    equation_b = sympy.Eq(sympy.Symbol('σ_b(ω)'), sigma_b)

    # Print the beautifully formatted equation
    sympy.pprint(equation_b, use_unicode=True)
    print("\nWhere the terms are defined as before, with ω_J being the new exciton transition frequency.")

if __name__ == '__main__':
    solve_absorption_cross_section()