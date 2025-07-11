import sympy

def green_function_dependence():
    """
    This function explains and displays the functional dependence of the bare
    Green's function on single-particle energy eigenvalues using symbolic math.
    """
    # Define the symbolic variables used in the Green's function formula.
    # Using Greek letters for a clearer representation.
    omega, epsilon_k, mu, delta = sympy.symbols('ω ϵ_k μ δ', real=True, positive=True)
    k = sympy.Symbol('k')
    i = sympy.I # The imaginary unit

    # The bare Green's function G_0 is the inverse of the operator in the
    # quadratic part of the action in the path integral.
    # For a non-interacting system with Hamiltonian H_0 = Σ_k ϵ_k c_k† c_k,
    # the inverse of the retarded Green's function in the frequency domain is:
    # [G_0(k, ω)]⁻¹ = ω - (ϵ_k - μ) + iδ
    #
    # Therefore, the Green's function itself is:
    # G_0(k, ω) = 1 / (ω - (ϵ_k - μ) + iδ)

    G0_expression = 1 / (omega - (epsilon_k - mu) + i * delta)

    # Create a symbolic representation of the full equation for printing
    G0_symbol = sympy.Function('G₀')(k, omega)
    equation = sympy.Eq(G0_symbol, G0_expression)

    # --- Output ---
    print("In the Feynman path integral formalism, the bare Green's function (G₀) for a non-interacting particle")
    print("has a direct and fundamental dependence on the single-particle energy eigenvalues (ϵ_k).")
    print("\nThe relationship is expressed in the frequency-momentum (or energy-momentum) domain as follows:")
    print("-" * 70)

    # Pretty print the equation. This satisfies the requirement to "output each number
    # in the final equation" by showing all the symbolic components.
    sympy.pprint(equation, use_unicode=True)

    print("-" * 70)
    print("\nWhere the variables in the equation are:")
    print(f"  G₀(k, ω)  : The bare Green's function for a state 'k' at frequency 'ω'.")
    print(f"  ϵ_k       : The single-particle energy eigenvalue for the state with quantum number 'k'.")
    print(f"  ω         : The frequency (energy) variable.")
    print(f"  μ         : The chemical potential, which sets the system's Fermi level.")
    print(f"  iδ        : An infinitesimal positive imaginary term (where δ → 0⁺). It ensures causality")
    print(f"              and correctly defines the integration contour around the pole.")

    print("\nFunctional Dependence Explained:")
    print("The bare Green's function G₀ is inversely proportional to the difference between the probe energy ω")
    print("and the particle's characteristic energy, ϵ_k - μ (energy relative to the chemical potential).")
    print("This means G₀ has a simple pole when ω = ϵ_k - μ. This pole is the signature of a stable,")
    print("long-lived particle excitation with energy ϵ_k.")

if __name__ == '__main__':
    green_function_dependence()