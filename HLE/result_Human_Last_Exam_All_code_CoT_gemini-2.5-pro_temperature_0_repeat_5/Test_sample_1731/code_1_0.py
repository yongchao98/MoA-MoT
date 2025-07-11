import sympy as sp

def derive_bose_einstein_thermodynamics():
    """
    Derives and prints the equilibrium values for mean energy and entropy
    for a system of photons (Bose-Einstein statistics) using symbolic math.

    The derivation is based on maximizing the entropy of the system, a principle
    related to the large deviation theorems mentioned in the problem.
    """
    # Define the symbolic variables for the equations.
    # β (beta) = 1 / (k_B * T) is the inverse temperature, where k_B is the
    # Boltzmann constant and T is the temperature.
    # ε_i is the energy of the i-th level.
    # g_i is the degeneracy of the i-th level (i.e., the number of states at that energy).
    beta, epsilon_i, g_i, k_B = sp.symbols('β ε_i g_i k_B', real=True, positive=True)

    # --- Step 1: The Bose-Einstein Distribution ---
    # The equilibrium (most probable) average occupation number <n_i> for a single state
    # at energy ε_i is found by maximizing the entropy for a system of bosons
    # with a fixed total energy and non-conserved particle number (the case for photons).
    # This is the core result from applying the entropy maximization principle.
    # The average number of particles per state is denoted <n_i> = n_i / g_i.
    avg_occupation_number = 1 / (sp.exp(beta * epsilon_i) - 1)

    print("--- Equilibrium Values for the Bose Case of Light Quanta ---")
    print("\nBased on the maximization of entropy, the equilibrium average occupation")
    print("number for a single state at energy ε_i is given by the Bose-Einstein distribution:")
    print("\n<n_i> = n_i / g_i =")
    sp.pprint(avg_occupation_number, use_unicode=True)
    print("\nwhere:")
    print("  n_i: number of photons in energy level ε_i")
    print("  g_i: number of degenerate states in energy level ε_i")
    print("  β: 1 / (k_B * T), the inverse temperature")
    print("-" * 60)

    # --- Step 2: Equilibrium Mean Energy (E) ---
    # The total mean energy of the system is the sum of the energies of all
    # photons. This is found by summing the energy of each level (ε_i) multiplied
    # by the total number of photons in that level (n_i = g_i * <n_i>).
    # The total energy E is the sum of contributions E_i from each level.
    # E = Σ_i E_i = Σ_i (n_i * ε_i)
    
    energy_contribution_i = g_i * avg_occupation_number * epsilon_i
    
    print("\n1. Equilibrium Mean Energy (E)")
    print("\nThe total mean energy E is the sum over all energy levels i: E = Σ_i E_i")
    print("The contribution to the mean energy from level 'i' is E_i = g_i * <n_i> * ε_i.")
    print("\nFinal Equation for the energy contribution from level i (E_i):")
    sp.pprint(energy_contribution_i, use_unicode=True)
    print("-" * 60)

    # --- Step 3: Equilibrium Entropy (S) ---
    # The entropy is calculated from the number of microstates W using S = k_B * ln(W).
    # For bosons, after using Stirling's approximation, the entropy contribution from level i is:
    # S_i / k_B = g_i * [ (<n_i> + 1)ln(<n_i> + 1) - <n_i>ln(<n_i>) ]
    # Substituting the equilibrium occupation number <n_i> gives the final form.
    # The total entropy S is the sum of contributions S_i from each level.
    
    entropy_contribution_i_over_kB = g_i * (
        (avg_occupation_number + 1) * sp.log(avg_occupation_number + 1) -
        avg_occupation_number * sp.log(avg_occupation_number)
    )
    
    # The expression can be simplified to show its relation to energy and the partition function.
    # S_i / k_B = β*E_i - g_i*ln(1 - exp(-β*ε_i))
    simplified_entropy_i = sp.simplify(entropy_contribution_i_over_kB)

    print("\n2. Equilibrium Entropy (S)")
    print("\nThe total entropy S is the sum over all energy levels i: S = Σ_i S_i")
    print("The contribution to the entropy from level 'i' is derived from the combinatorial formula.")
    print("\nFinal Equation for the entropy contribution from level i (S_i / k_B):")
    sp.pprint(simplified_entropy_i, use_unicode=True)
    print("-" * 60)


if __name__ == '__main__':
    derive_bose_einstein_thermodynamics()