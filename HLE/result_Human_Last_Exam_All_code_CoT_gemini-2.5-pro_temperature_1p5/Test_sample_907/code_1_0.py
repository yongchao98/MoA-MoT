import sympy

def display_absorption_equations():
    """
    This function defines and prints the equations for the absorption cross-section
    for a chain of molecules under two different assumptions, using sympy for
    mathematical representation.
    """
    # Define the symbols used in the equations
    sigma_abs = sympy.Symbol('σ_abs')
    C = sympy.Symbol('C')  # Proportionality constant
    N = sympy.Symbol('N')  # Number of molecules in the chain
    mu0 = sympy.Symbol('μ₀')  # Transition dipole moment of a single molecule
    E_trans = sympy.Symbol('E_trans')  # Excitation energy of a single molecule
    hbar = sympy.Symbol('ħ') # Reduced Planck constant
    omega_L = sympy.Symbol('ω_L')  # Central frequency of the laser pulse
    tau = sympy.Symbol('τ')  # Duration of the laser pulse
    J = sympy.Symbol('J')  # Near-neighbor coupling energy

    # --- Case a) No interaction between molecules ---
    # The total absorption is the sum of N identical, non-interacting molecules.
    # The transition dipole moment squared for the system is N * |μ₀|².
    # The transition energy is just E_trans.
    exponent_a = -(((E_trans - hbar * omega_L) * tau / hbar)**2)
    equation_a = sympy.Eq(sigma_abs, C * N * sympy.Abs(mu0)**2 * sympy.exp(exponent_a))

    # --- Case b) Near-neighbor interaction ---
    # Interaction leads to delocalized exciton states. Only the k=0 state is optically active.
    # The energy of the k=0 exciton state is shifted by 2J.
    # The oscillator strength is coherently summed into this one state, giving a transition
    # dipole moment squared of N * |μ₀|².
    exponent_b = -(((E_trans + 2 * J - hbar * omega_L) * tau / hbar)**2)
    equation_b = sympy.Eq(sigma_abs, C * N * sympy.Abs(mu0)**2 * sympy.exp(exponent_b))

    # --- Print the final equations in a readable format ---
    print("The equation for the absorption cross-section (σ_abs) for a chain of molecules is given by:")
    print("\nSymbols represent:")
    print(f"  C: Proportionality constant")
    print(f"  N: Number of molecules")
    print(f"  μ₀: Transition dipole moment of a single molecule")
    print(f"  E_trans: Transition energy of a single molecule")
    print(f"  ħ: Reduced Planck constant")
    print(f"  ω_L: Laser pulse central frequency")
    print(f"  τ: Laser pulse duration")
    print(f"  J: Near-neighbor coupling energy")

    print("\n----------------------------------------------------------")
    print("a) The interaction between molecules can be neglected:")
    print("----------------------------------------------------------")
    sympy.pprint(equation_a, use_unicode=True)
    print("\nHere, the absorption peak is centered at the single-molecule transition energy E_trans.")

    print("\n----------------------------------------------------------")
    print("b) The interaction between near-neighbors should be considered:")
    print("----------------------------------------------------------")
    sympy.pprint(equation_b, use_unicode=True)
    print("\nHere, the absorption peak is shifted to E_trans + 2J due to exciton formation.")


if __name__ == '__main__':
    display_absorption_equations()