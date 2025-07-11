import sympy

def solve_bose_equilibrium():
    """
    Derives and prints the equilibrium values for mean energy and entropy
    for a gas of light quanta (Bose case) using the principle of maximum entropy.
    """

    # --- Introduction ---
    print("Derivation of Equilibrium Properties for a Photon Gas (Bose case)")
    print("------------------------------------------------------------------")
    print("We use the principle of maximum entropy, a concept at the heart of")
    print("Boltzmann-Sanov and Chernoff large deviation theorems.")
    print("The equilibrium state is the one that maximizes the entropy S for a given mean energy U.")
    print("\nFor bosons (like photons) in an energy level 'i' with energy epsilon_i,")
    print("the number of microstates Omega_i for N_i particles in g_i states is:")
    print("Omega_i = (N_i + g_i - 1)! / (N_i! * (g_i - 1)!)")
    print("The entropy is S_i = k_B * ln(Omega_i). We use Stirling's approximation for ln(n!).")
    print("(For simplicity, we set Boltzmann's constant k_B = 1)\n")


    # --- Step 1: Define Symbols and Entropy ---
    # N_i: number of photons in level i
    # g_i: number of states (degeneracy) in level i
    # epsilon_i: energy of level i
    # beta: Lagrange multiplier, physically 1/(k_B*T)
    # <n_i>: mean occupation number per state in level i
    N_i, g_i, epsilon_i, beta = sympy.symbols('N_i g_i epsilon_i beta', positive=True, real=True)
    n_i = sympy.symbols('<n_i>', positive=True, real=True)

    # Entropy S_i using Stirling's approximation: S ≈ (N+g)ln(N+g) - Nln(N) - gln(g)
    S_i_approx = (N_i + g_i) * sympy.log(N_i + g_i) - N_i * sympy.log(N_i) - g_i * sympy.log(g_i)

    # --- Step 2: Maximize Entropy with Lagrange Multiplier ---
    # We maximize S_i with the constraint of fixed mean energy U_i = N_i * epsilon_i.
    # We form the Lagrangian L = S_i - beta * U_i and find where dL/dN_i = 0.
    # This leads to the condition: dS_i/dN_i = beta * epsilon_i
    dS_dNi = sympy.diff(S_i_approx, N_i)
    # The derivative dS_dNi simplifies to ln((N_i + g_i) / N_i)
    # So we solve: ln((N_i + g_i) / N_i) = beta * epsilon_i
    # Exponentiating both sides: (N_i + g_i) / N_i = exp(beta * epsilon_i)
    # 1 + g_i / N_i = exp(beta * epsilon_i)
    # g_i / N_i = exp(beta * epsilon_i) - 1
    # N_i / g_i = 1 / (exp(beta * epsilon_i) - 1)
    # This is the mean occupation number per state, <n_i>
    n_i_equilibrium = 1 / (sympy.exp(beta * epsilon_i) - 1)

    print("--- 1. Equilibrium Occupation Number ---")
    print("By maximizing entropy, we find the equilibrium mean occupation number <n_i> for a state with energy epsilon_i:")
    print("\n<n_i> =")
    sympy.pprint(n_i_equilibrium)
    print("\nThis is the Bose-Einstein distribution for photons (chemical potential = 0).")
    print("\nWhere:")
    print("  epsilon_i : Energy of a single quantum state.")
    print("  beta      : Inverse temperature (1 / (k_B * T)).\n")


    # --- Step 3: Define and Print Mean Energy Equation ---
    # The total mean energy U is the sum over all states of the mean energy per state.
    # Mean energy per state is <n_i> * epsilon_i
    mean_energy_per_state = n_i * epsilon_i
    print("--- 2. Equilibrium Mean Energy (U) ---")
    print("The total mean energy U of the system is the sum of the energies of all occupied states.")
    print("The contribution from each state with energy epsilon_i is <n_i> * epsilon_i.")
    print("So, the total mean energy is given by the integral (or sum) over all states:")
    print("\nU = Sum_i [ <n_i> * epsilon_i ] = Sum_i [")
    sympy.pprint(n_i_equilibrium * epsilon_i)
    print("]")
    print("(The sum is over all energy levels i, often replaced by an integral over the density of states).\n")

    # --- Step 4: Define and Print Entropy Equation ---
    # The equilibrium entropy per state, s_i = S_i/g_i, is given by the formula:
    # s_i = (<n_i>+1)ln(<n_i>+1) - <n_i>ln(<n_i>)
    entropy_per_state = (n_i + 1) * sympy.log(n_i + 1) - n_i * sympy.log(n_i)
    print("--- 3. Equilibrium Entropy (S) ---")
    print("The total entropy S is the sum of the entropy contributions from each energy level.")
    print("The entropy per state is expressed in terms of the mean occupation number <n_i>:")
    print("\nS = Sum_i [ s_i ] = Sum_i [")
    sympy.pprint(entropy_per_state)
    print("]")
    print("\nTo find the value of S, one first calculates <n_i> for each level using the Bose-Einstein formula,")
    print("and then substitutes that value into the expression for s_i above.")


if __name__ == '__main__':
    solve_bose_equilibrium()
    # The final answer is the set of derived equations for U and S.
    # The problem asks for values, which in this theoretical context are the functional forms.
    # We select the entropy formula as the final representative answer.
    final_answer = "S = Sum_i [ (<n_i>+1)*ln(<n_i>+1) - <n_i>*ln(<n_i>) ] where <n_i> = 1/(exp(beta*epsilon_i)-1)"
    # print(f"<<<{final_answer}>>>")