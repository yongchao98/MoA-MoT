import sympy as sp

def solve_bose_light_quanta():
    """
    Derives the equilibrium values of mean energy (U) and entropy (S)
    for a gas of light quanta (photons) using principles of statistical mechanics.
    """
    # --- Step 1: Define physical constants and variables symbolically ---
    # We use sympy for symbolic mathematics.
    epsilon, T, V = sp.symbols('epsilon T V', positive=True, real=True)
    k, h, c = sp.symbols('k h c', positive=True, real=True) # Boltzmann, Planck, speed of light
    
    # The Lagrange multiplier beta is related to temperature T
    beta = 1 / (k * T)

    print("Derivation of Equilibrium Values for Light Quanta (Bose Case)")
    print("="*60)

    # --- Step 2: State the Bose-Einstein Distribution ---
    # The equilibrium occupation number n(epsilon) for a state with energy epsilon is found
    # by maximizing the system's entropy for a fixed mean energy. The result
    # is the Bose-Einstein distribution. For photons, the chemical potential is zero.
    n_epsilon = 1 / (sp.exp(beta * epsilon) - 1)
    
    print("\n1. The Bose-Einstein distribution n(ε) for photons is:")
    sp.pprint(sp.Eq(sp.Symbol('n(ε)'), n_epsilon))

    # --- Step 3: Define the Density of States ---
    # For photons (relativistic bosons) in a 3D box of volume V, the number of
    # available states per unit energy is g(epsilon).
    g_epsilon = (8 * sp.pi * V / (h**3 * c**3)) * epsilon**2
    
    print("\n2. The density of states g(ε) for photons in a volume V is:")
    sp.pprint(sp.Eq(sp.Symbol('g(ε)'), g_epsilon))

    # --- Step 4: Calculate the Mean Energy (U) ---
    # U is the integral of (energy * occupation_number * density_of_states) over all energies.
    # The integral to solve is ∫ ε * n(ε) * g(ε) dε from 0 to ∞.
    integrand_U = epsilon * n_epsilon * g_epsilon
    
    # We can ask sympy to perform the definite integral from 0 to infinity.
    # sympy knows the value of the standard Bose-Einstein integral.
    U = sp.integrate(integrand_U, (epsilon, 0, sp.oo))
    
    print("\n3. The equilibrium mean energy U is calculated by integrating ε*n(ε)*g(ε):")
    U_integral_form = sp.Integral(integrand_U, (epsilon, 0, sp.oo))
    sp.pprint(sp.Eq(sp.Symbol('U'), U_integral_form))
    
    print("\n   Evaluating the integral gives the Stefan-Boltzmann Law for total energy:")
    sp.pprint(sp.Eq(sp.Symbol('U'), U))

    # --- Step 5: Calculate the Entropy (S) ---
    # For a photon gas (an ultra-relativistic Bose gas), there is a simple thermodynamic
    # relationship between entropy S, mean energy U, and temperature T.
    # This can be derived from the expression for the Helmholtz free energy.
    # The relationship is S = (4/3) * U / T
    S = (4 * U) / (3 * T)
    
    print("\n4. The equilibrium entropy S for a photon gas is related to U and T:")
    sp.pprint(sp.Eq(sp.Symbol('S'), S.simplify()))

    # --- Step 6: Final Equation Output ---
    # The question asks to print the final equation with each part separated.
    print("\n5. The final relationship between equilibrium entropy and mean energy is:")
    
    # The equation is S = (4/3) * U / (k*T). Let's print it piece by piece.
    final_equation_parts = ["S", "=", "4", "/", "3", "*", "U", "/", "(", "k", "*", "T", ")"]
    print(" ".join(final_equation_parts))

if __name__ == '__main__':
    solve_bose_light_quanta()