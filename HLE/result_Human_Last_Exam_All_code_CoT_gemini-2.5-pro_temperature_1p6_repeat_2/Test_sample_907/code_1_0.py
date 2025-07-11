import sympy as sp

def display_absorption_equations():
    """
    This function uses the sympy library to symbolically construct and print
    the equations for the absorption cross-section of a molecular chain.
    """

    # --- Define Symbolic Variables ---
    # Using unicode characters for a clear mathematical representation.
    C_const = sp.Symbol('C', real=True, positive=True)
    omega = sp.Symbol('ω', real=True, positive=True)
    tau = sp.Symbol('τ', real=True, positive=True)
    omega_eg = sp.Symbol('ω_eg', real=True, positive=True)
    mu_eg_sq = sp.Symbol('|μ_eg|²', real=True, positive=True)
    N = sp.Symbol('N', integer=True, positive=True)
    J = sp.Symbol('J', real=True)
    hbar = sp.Symbol('ħ', real=True, positive=True)
    sigma_a = sp.Function('σ_a')(omega)
    sigma_b = sp.Function('σ_b')(omega)

    print("This program provides the equations for the absorption cross-section (σ) of a molecular chain.")
    print("The system is excited by an ultrashort Gaussian laser pulse of duration τ.")
    print("\nThe general form of the equation is: σ(ω) ∝ (Transition Strength) ⋅ ω ⋅ (Lineshape)")
    print("-" * 70)

    # --- Case a) No interaction between molecules ---
    print("\na) The interaction between molecules can be neglected.\n")

    # Define each component of the final equation for clarity
    strength_a = N * mu_eg_sq
    frequency_a = omega_eg
    lineshape_a = sp.exp(-(omega - frequency_a)**2 * tau**2)

    # Combine components into the final equation
    equation_a = C_const * strength_a * omega * lineshape_a

    print("The components of the equation are:")
    sp.pprint(sp.Eq(sp.Symbol('Transition Strength'), strength_a), use_unicode=True)
    sp.pprint(sp.Eq(sp.Symbol('Transition Frequency'), frequency_a), use_unicode=True)
    sp.pprint(sp.Eq(sp.Symbol('Gaussian Lineshape'), lineshape_a), use_unicode=True)
    print("\nThe final equation for the absorption cross-section σ_a(ω) is:")
    sp.pprint(sp.Eq(sigma_a, equation_a), use_unicode=True)
    print("\nWhere C is a proportionality constant, N is the number of molecules, |μ_eg|² is the monomer's\n"
          "squared transition dipole moment, ω is the photon frequency, ω_eg is the monomer's transition\n"
          "frequency, and τ is the laser pulse duration.")

    print("-" * 70)

    # --- Case b) Nearest-neighbor interaction ---
    print("\nb) The interaction between near-neighbors should be considered.\n")

    # Define each component of the final equation for clarity
    strength_b = N * mu_eg_sq  # Coherently enhanced strength
    frequency_b = omega_eg + (2 * J / hbar) # Shifted frequency due to exciton coupling
    lineshape_b = sp.exp(-(omega - frequency_b)**2 * tau**2)

    # Combine components into the final equation
    equation_b = C_const * strength_b * omega * lineshape_b

    print("The components of the equation are:")
    sp.pprint(sp.Eq(sp.Symbol('Transition Strength'), strength_b), use_unicode=True)
    sp.pprint(sp.Eq(sp.Symbol('Exciton Transition Frequency'), frequency_b), use_unicode=True)
    sp.pprint(sp.Eq(sp.Symbol('Gaussian Lineshape'), lineshape_b), use_unicode=True)
    print("\nThe final equation for the absorption cross-section σ_b(ω) is:")
    sp.pprint(sp.Eq(sigma_b, equation_b), use_unicode=True)
    print("\nWhere J is the nearest-neighbor interaction energy and ħ is the reduced Planck constant.\n"
          "The interaction shifts the absorption peak by 2J/ħ and coherently enhances its strength.")

if __name__ == '__main__':
    display_absorption_equations()