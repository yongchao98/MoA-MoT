import sympy

def solve_bose_einstein_equilibrium():
    """
    This script derives the equilibrium values for mean energy and entropy
    for a gas of light quanta (photons) using symbolic mathematics.

    The physical principle, supported by large deviation theorems (like Boltzmann-Sanov),
    is that the equilibrium state of an isolated system is the one that maximizes
    entropy subject to physical constraints.

    For a photon gas (Bose statistics):
    1. Entropy (S): S = k * ln(W), where W is the number of microstates.
       Using Stirling's approximation, S/k becomes a sum over energy levels 'i'.
    2. Constraint: The total energy E = Σ (n_i * ε_i) is constant.
    3. Non-Constraint: The number of photons N = Σ n_i is NOT conserved.

    We use the method of Lagrange multipliers to maximize S subject to the energy constraint.
    """
    print("--- Derivation of the Equilibrium Distribution ---")
    print("We maximize the function L = S/k - β * E, where β is a Lagrange multiplier.")
    print("This is equivalent to solving ∂(S/k - βE)/∂n_i = 0 for each energy level i.")

    # Define symbolic variables using sympy
    # Using English names for variables for broader compatibility
    n, g, epsilon, beta = sympy.symbols('n_i g_i epsilon_i beta', real=True, positive=True)

    # Entropy term for a single energy level 'i' (S_i / k)
    # Using the approximation for large numbers: ln(W_i) ≈ (n+g)ln(n+g) - n*ln(n) - g*ln(g)
    log_W_i = (n + g) * sympy.log(n + g) - n * sympy.log(n) - g * sympy.log(g)

    # Lagrangian for a single level 'i'
    Lagrangian_i = log_W_i - beta * n * epsilon
    
    # Use pretty symbols for printing
    n_sym, g_sym, e_sym, b_sym = sympy.symbols('n_i g_i ε_i β')
    lag_print = ((n_sym + g_sym) * sympy.log(n_sym + g_sym) - n_sym * sympy.log(n_sym) - g_sym * sympy.log(g_sym)) - b_sym*n_sym*e_sym


    print("\nThe Lagrangian term for a single energy level 'i' is:")
    print(f"L_i = S_i/k - β*n_i*ε_i = {sympy.pretty(lag_print, use_unicode=True)}")

    # Differentiate with respect to n
    dL_dn = sympy.diff(Lagrangian_i, n)
    dL_dn_print = sympy.diff(lag_print, n_sym)
    
    print(f"\nDifferentiating with respect to n_i and setting to 0:")
    print(f"∂L_i/∂n_i = {sympy.pretty(dL_dn_print, use_unicode=True)} = 0")

    # Solve the equation dL/dn = 0 for n
    equilibrium_eq = sympy.Eq(dL_dn, 0)
    solution = sympy.solve(equilibrium_eq, n)
    n_equilibrium = solution[0]

    print("\nSolving for n_i gives the equilibrium occupation number (Bose-Einstein distribution):")
    # Using f-string to insert numbers like '1' explicitly
    final_dist_str = f"n_i = g_i / (exp(β*ε_i) - 1)"
    print(final_dist_str)
    print("where g_i is the number of states at energy level ε_i, and the Lagrange multiplier β is related to temperature T by β = 1/(k*T).")
    print("The average occupation per state is <n_i> = n_i/g_i.")


    print("\n\n--- Equilibrium Mean Energy (E_eq) ---")
    print("The mean energy of the system is found by summing the energy of all occupied states using the equilibrium distribution:")
    # Again using f-string to ensure numbers are present
    energy_eq_str = "E_eq = Σ_i [ n_i * ε_i ] = Σ_i [ g_i * ε_i / (exp(β*ε_i) - 1) ]"
    print(energy_eq_str)


    print("\n\n--- Equilibrium Entropy (S_eq) ---")
    print("The equilibrium entropy is found by substituting the equilibrium distribution n_i back into the entropy formula.")
    print("This yields the following expression, expressed in terms of the average occupation number <n_i> = n_i/g_i:")
    entropy_eq_str = "S_eq = k * Σ_i g_i * [ (1 + <n_i>) * ln(1 + <n_i>) - <n_i> * ln(<n_i>) ]"
    print(entropy_eq_str)

    print("\nAn alternative (and equivalent) expression is derived from the thermodynamic relation S = (E - A)/T, where A is the Helmholtz free energy:")
    entropy_alt_str = "S_eq = k * [ β*E_eq - Σ_i g_i * ln(1 - exp(-β*ε_i)) ]"
    print(entropy_alt_str)


if __name__ == '__main__':
    solve_bose_einstein_equilibrium()