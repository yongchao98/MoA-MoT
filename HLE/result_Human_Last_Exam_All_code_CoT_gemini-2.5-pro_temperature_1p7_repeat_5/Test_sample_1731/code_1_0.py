import sympy

def solve_bose_einstein_equilibrium():
    """
    Derives the equilibrium mean energy and entropy for a Bose gas (photons)
    using the principle of maximum entropy.
    """

    # --- Introduction ---
    print("Derivation of Equilibrium Energy and Entropy for a Bose Gas (Photons)")
    print("="*70)
    print("We will use the principle of maximum entropy, which is the physical application of")
    print("large deviation theory, to find the most probable (equilibrium) state.")
    print("\nLet's consider a single energy level 'i' with:")
    print(" - n_i: number of particles (bosons)")
    print(" - g_i: number of degenerate states (degeneracy)")
    print(" - epsilon_i: energy of the level")
    print(" - k_B: Boltzmann constant")
    print(" - beta: Lagrange multiplier, which will be identified with 1/(k_B*T)")
    print("-" * 70, "\n")

    # --- Symbolic Setup ---
    # Define symbols for a generic energy level i
    n_i, g_i, epsilon_i = sympy.symbols('n_i g_i epsilon_i', real=True, positive=True)
    beta, k_B = sympy.symbols('beta k_B', real=True, positive=True)

    # --- Step 1: Entropy and Stirling's Approximation ---
    print("Step 1: Define Entropy from Combinatorics")
    print("The number of ways W_i to place n_i indistinguishable bosons in g_i states is:")
    print("W_i = (n_i + g_i - 1)! / (n_i! * (g_i - 1)!)")
    print("\nThe entropy S is S = k_B * ln(W), where W = Product_i(W_i).")
    print("We maximize S_i/k_B = ln(W_i). For large n_i and g_i, we use Stirling's approximation:")
    print("ln(x!) ≈ x*ln(x) - x")
    
    # Entropy for level i (divided by k_B) after applying Stirling's approximation
    log_W_i = (n_i + g_i) * sympy.log(n_i + g_i) - n_i * sympy.log(n_i) - g_i * sympy.log(g_i)
    print("\nApplying this, the entropy contribution from level i is approximately:")
    sympy.pprint(log_W_i, use_unicode=True)
    print("-" * 70, "\n")
    
    # --- Step 2: Maximization using Lagrange Multipliers ---
    print("Step 2: Maximize Entropy with Energy Constraint")
    print("We maximize the total entropy Sum(S_i) subject to constant total energy Sum(n_i * epsilon_i) = E.")
    print("We form the Lagrangian L = S/k_B - beta * E = Sum_i (ln(W_i) - beta * n_i * epsilon_i).")
    print("To find the maximum, we set the derivative dL/dn_i = 0.")
    
    # Define the Lagrangian for level i
    L_i = log_W_i - beta * n_i * epsilon_i
    
    # Differentiate with respect to n_i
    deriv_L_i = sympy.diff(L_i, n_i)
    
    print("\nThe derivative of the Lagrangian with respect to n_i is:")
    sympy.pprint(deriv_L_i, use_unicode=True)

    # Solve for n_i
    equilibrium_n_i_expr = sympy.solve(deriv_L_i, n_i)[0]
    
    print("\nSolving dL/dn_i = 0 for n_i gives the equilibrium number of particles in level i:")
    sympy.pprint(sympy.Eq(n_i, equilibrium_n_i_expr), use_unicode=True)
    print("-" * 70, "\n")
    
    # --- Step 3: Equilibrium Distribution and Energy ---
    print("Step 3: Find Equilibrium Occupation Number and Mean Energy")
    print("The average occupation number <n_i> is the number of particles per state, n_i/g_i.")
    avg_occupation = equilibrium_n_i_expr / g_i
    
    n_i_avg = sympy.Symbol('<n_i>')
    print("\nThe equilibrium average occupation number <n_i> is therefore:")
    final_occupation_eq = sympy.Eq(n_i_avg, avg_occupation)
    sympy.pprint(final_occupation_eq, use_unicode=True)
    print("\nThis is the Bose-Einstein distribution for photons (chemical potential = 0).")
    
    print("\nThe equilibrium mean energy E is the sum of energies over all levels:")
    E_symbol = sympy.Symbol('E')
    epsilon_i_symbol = sympy.Symbol('epsilon_i')
    g_i_symbol = sympy.Symbol('g_i')
    # Using the symbol for <n_i> to present the formula
    E_eq = sympy.Sum(g_i_symbol * epsilon_i_symbol * n_i_avg, (sympy.Symbol('i'), 1, sympy.Symbol('M')))
    
    print("\nEQUILIBRIUM MEAN ENERGY:")
    sympy.pprint(sympy.Eq(E_symbol, E_eq), use_unicode=True)
    print("where <n_i> is given by the distribution above.")
    print("-" * 70, "\n")

    # --- Step 4: Equilibrium Entropy ---
    print("Step 4: Find Equilibrium Entropy")
    print("Substituting the equilibrium n_i back into the entropy formula S = k_B * Sum(ln(W_i))")
    print("and simplifying gives the final expression for the equilibrium entropy.")

    # The simplified form of the entropy for a Bose gas
    S_symbol = sympy.Symbol('S')
    entropy_term_i = g_i * ((n_i_avg + 1) * sympy.log(n_i_avg + 1) - n_i_avg * sympy.log(n_i_avg))
    S_eq = k_B * sympy.Sum(entropy_term_i.subs(g_i, g_i_symbol), (sympy.Symbol('i'), 1, sympy.Symbol('M')))
    
    print("\nEQUILIBRIUM ENTROPY:")
    sympy.pprint(sympy.Eq(S_symbol, S_eq), use_unicode=True)
    print("where <n_i> is the Bose-Einstein distribution function derived above.")
    print("="*70)

# Run the derivation
solve_bose_einstein_equilibrium()