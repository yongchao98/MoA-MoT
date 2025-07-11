import sympy as sp

def derive_bose_einstein_equilibrium():
    """
    Derives the equilibrium mean energy and entropy for a Bose gas of photons
    by maximizing the Boltzmann entropy.
    """
    # 1. Define the necessary symbolic variables.
    # n_i: occupation number of energy level i
    # g_i: degeneracy (number of states) of energy level i
    # e_i: energy of level i
    # k: Boltzmann constant
    # T: Temperature
    ni, gi, ei = sp.symbols('n_i g_i e_i', real=True, positive=True)
    k, T = sp.symbols('k T', real=True, positive=True)

    # The Lagrange multiplier beta is related to temperature by beta = 1/(k*T).
    beta = sp.Symbol('beta', real=True, positive=True)

    print("Step 1: Define the entropy for a single energy level `i` for a Bose gas.")
    print("Using Stirling's approximation on the Bose-Einstein counting formula, the entropy contribution S_i is:")
    # S_i/k = (n_i + g_i)ln(n_i + g_i) - n_i*ln(n_i) - g_i*ln(g_i)
    # This is the quantity we need to maximize.
    S_over_k_i = (ni + gi) * sp.log(ni + gi) - ni * sp.log(ni) - gi * sp.log(gi)
    print(f"S_i / k = {S_over_k_i}\n")

    print("Step 2: Maximize entropy subject to energy conservation using Lagrange multipliers.")
    print("For photons, particle number is not conserved, so we only have one constraint: sum(n_i*e_i) = E.")
    print("We maximize L = S - beta * E. This means solving d(S_i/k - beta*n_i*e_i)/d(n_i) = 0 for n_i.")
    
    # The derivative of the Lagrangian term for level i with respect to n_i
    lagrangian_derivative = sp.diff(S_over_k_i, ni) - beta * ei
    
    print(f"Setting the derivative to zero: {sp.diff(S_over_k_i, ni)} - beta * e_i = 0\n")

    # 2. Solve for the equilibrium occupation number n_i.
    equilibrium_eq = sp.Eq(sp.diff(S_over_k_i, ni), beta * ei)
    ni_eq_sol = sp.solve(equilibrium_eq, ni)
    # The solution is gi / (exp(beta*ei) - 1)
    ni_eq = ni_eq_sol[0]
    
    # Substitute beta = 1/(k*T) for the final expression
    ni_eq_T = ni_eq.subs(beta, 1/(k*T))

    print("Step 3: Solve for the equilibrium occupation number n_i (Bose-Einstein distribution).")
    print("The equilibrium number of particles in a state with energy e_i is:")
    final_ni_eq = sp.Eq(ni, ni_eq_T)
    print(f"{sp.pretty(final_ni_eq, use_unicode=False)}\n")

    # 3. Calculate the equilibrium mean energy for level i.
    Ei_eq = ni_eq_T * ei

    print("Step 4: Calculate the equilibrium mean energy for level i.")
    print("The mean energy of level `i`, E_i, is the occupation number n_i times the energy per particle e_i.")
    final_Ei_eq = sp.Eq(sp.Symbol('E_i'), Ei_eq)
    print(f"{sp.pretty(final_Ei_eq, use_unicode=False)}\n")

    # 4. Calculate the equilibrium entropy for level i.
    # The entropy S_i can be expressed thermodynamically based on the derived n_i.
    # A standard and clean form is: S_i = E_i/T - k*g_i*ln(1 - exp(-e_i/kT))
    # Let's construct this expression.
    Si_eq = Ei_eq/T - k*gi*sp.log(1 - sp.exp(-ei/(k*T)))
    
    print("Step 5: Calculate the equilibrium entropy for level i.")
    print("The entropy contribution from level `i`, S_i, is found by substituting n_i back into the entropy formula.")
    print("A simplified and standard result from thermodynamics is:")
    final_Si_eq = sp.Eq(sp.Symbol('S_i'), Si_eq)
    print(f"{sp.pretty(final_Si_eq, use_unicode=False)}\n")
    
    print("---")
    print("Summary of Equilibrium Values:")
    print("The total energy E and entropy S of the photon gas are obtained by summing E_i and S_i over all possible states.")
    print(f"Equilibrium Mean Energy per state: {final_Ei_eq}")
    print(f"Equilibrium Entropy per state: {final_Si_eq}")


if __name__ == '__main__':
    derive_bose_einstein_equilibrium()