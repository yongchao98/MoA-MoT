import math

def explain_chemical_potential_limit():
    """
    Explains the fundamental limit on the chemical potential for bosons in a Bose-Einstein condensate.
    """
    print("1. The Bose-Einstein distribution gives the occupation number 'n_i' of a state with energy 'epsilon_i':")
    print("   n_i = 1 / (exp((epsilon_i - mu) / (k_B * T)) - 1)")
    print("-" * 50)

    print("2. For 'n_i' to be a positive number, the denominator must be positive.")
    print("   This means: exp((epsilon_i - mu) / (k_B * T)) > 1")
    print("   Which simplifies to the condition: mu < epsilon_i")
    print("-" * 50)

    print("3. This condition must hold for ALL energy states. The strictest limit is set by the lowest energy state, the ground state 'epsilon_0'.")
    print("   Therefore, the fundamental limit is: mu < epsilon_0")
    print("-" * 50)

    print("4. During Bose-Einstein condensation (at T <= T_c), a macroscopic number of bosons occupies the ground state.")
    print("   This forces the chemical potential 'mu' to become equal to the ground state energy 'epsilon_0'.")
    print("-" * 50)

    print("5. The chemical potential of a non-interacting Bose gas at absolute zero (T=0) is also equal to the ground state energy, as all particles are in that state.")
    print("   So, mu(T=0) = epsilon_0.")
    print("-" * 50)

    print("Final Conclusion: The chemical potential in a Bose-Einstein Condensate is equal to the ground state energy, which is also the chemical potential of the gas at T=0.")
    mu_bec = "mu_BEC"
    epsilon_0 = "epsilon_0"
    mu_T0 = "mu(T=0)"
    print(f"The resulting equation is:")
    print(f"{mu_bec} = {epsilon_0}")
    print(f"And since {epsilon_0} = {mu_T0}, the final relationship is:")
    print(f"{mu_bec} = {mu_T0}")

explain_chemical_potential_limit()