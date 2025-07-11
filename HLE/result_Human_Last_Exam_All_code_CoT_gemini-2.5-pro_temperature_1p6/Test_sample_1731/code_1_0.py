import sympy

def solve_photon_gas_equilibrium():
    """
    Derives and prints the equilibrium values for mean energy and entropy of a photon gas
    using principles of statistical mechanics, as requested.
    """
    # Define physical constants and variables as symbols for symbolic manipulation
    V, T, k_B, h, c = sympy.symbols('V T k_B h c', positive=True, real=True)
    e, x = sympy.symbols('e x', positive=True, real=True)

    # In statistical mechanics, beta is a convenient shorthand for 1/(k_B*T)
    beta = 1 / (k_B * T)

    # --- Part 1: Derivation of Equilibrium Mean Energy <E> ---
    print("Derivation of Equilibrium Mean Energy <E>")
    print("=" * 45)
    print("The equilibrium state maximizes entropy, which corresponds to the Bose-Einstein distribution.")
    
    # The density of states g(e) for photons in a 3D cavity of volume V
    g_e = (8 * sympy.pi * V / (h**3 * c**3)) * e**2
    print("\n1. The density of available states g(e) for a photon with energy 'e' is:")
    sympy.pretty_print(sympy.Eq(sympy.Symbol('g(e)'), g_e))

    # The Bose-Einstein distribution n(e) for photons (which have zero chemical potential)
    n_e = 1 / (sympy.exp(e * beta) - 1)
    print("\n2. The average number of photons n(e) in a state with energy 'e' is:")
    sympy.pretty_print(sympy.Eq(sympy.Symbol('n(e)'), n_e))

    # The mean energy <E> is the integral of energy * occupation_number * density_of_states
    integrand_E = e * n_e * g_e
    integral_expr_E = sympy.Integral(integrand_E, (e, 0, sympy.oo))
    print("\n3. The mean energy <E> is the integral over all energies:")
    sympy.pretty_print(sympy.Eq(sympy.Symbol('<E>'), integral_expr_E))

    # To solve the integral, we perform a substitution x = beta*e. The integral part
    # becomes integral[x^3 / (exp(x) - 1) dx], which is a standard definite integral
    # related to the Riemann zeta function and evaluates to pi**4 / 15.
    
    # Final assembly of the expression for <E>
    C_E_part = (8 * sympy.pi * V / (h**3 * c**3))
    integral_E_val = sympy.pi**4 / 15
    mean_E = C_E_part * (k_B*T)**4 * integral_E_val
    mean_E_simplified = sympy.simplify(mean_E)

    print("\n4. After solving the integral, the final equation for mean energy is:")
    sympy.pretty_print(sympy.Eq(sympy.Symbol('<E>'), mean_E_simplified))

    print("\nBreaking down the numbers in the final equation for energy:")
    print("  <E> = (8 * pi^5 * V * k_B^4 * T^4) / (15 * h^3 * c^3)")
    print("  The specific numbers and powers are:")
    print("  - Numerator: 8, 5 (power of pi), 4 (power of k_B and T)")
    print("  - Denominator: 15, 3 (power of h and c)")
    

    # --- Part 2: Derivation of Equilibrium Entropy S ---
    print("\n\nDerivation of Equilibrium Entropy S")
    print("=" * 45)
    
    # For a photon gas (an ultra-relativistic gas), S = (4/3) * <E> / T.
    print("1. For a photon gas, entropy S is related to energy <E> by the simple thermodynamic formula:")
    print("   S = (4/3) * <E> / T")

    # Substitute the derived expression for <E> to find S
    entropy_S = (sympy.S(4)/3) * mean_E_simplified / T
    entropy_S_simplified = sympy.simplify(entropy_S)

    print("\n2. Substituting the result for <E>, we find the equilibrium entropy:")
    sympy.pretty_print(sympy.Eq(sympy.Symbol('S'), entropy_S_simplified))

    print("\nBreaking down the numbers in the final equation for entropy:")
    print("  S = (32 * pi^5 * V * k_B^4 * T^3) / (45 * h^3 * c^3)")
    print("  The specific numbers and powers are:")
    print("  - Numerator: 32 (from 4*8), 5 (power of pi), 4 (power of k_B)")
    print("  - Denominator: 45 (from 3*15), 3 (power of h, c and T)")


if __name__ == '__main__':
    solve_photon_gas_equilibrium()