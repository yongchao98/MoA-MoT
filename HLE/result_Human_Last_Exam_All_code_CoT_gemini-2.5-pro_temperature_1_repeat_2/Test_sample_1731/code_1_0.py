import sympy as sp

def find_equilibrium_values():
    """
    Derives and prints the equilibrium values of entropy and mean energy
    for the Bose case of light quanta (photons) based on statistical mechanics,
    which is underpinned by large deviation principles.
    """

    # --- Introduction and Connection to Large Deviation Theory ---
    print("--- Equilibrium Properties of a Photon Gas (Bose Case) ---")
    print("\nPrinciple: The equilibrium state of a physical system corresponds to the macrostate with the maximum entropy (S)")
    print("for a given set of constraints, such as a fixed mean energy (E). This is a core concept in statistical")
    print("mechanics, directly supported by large deviation theorems (e.g., Boltzmann-Sanov, Cramer-Chernoff).")
    print("These theorems mathematically show that the probability of observing a system in a state significantly")
    print("different from this maximum-entropy equilibrium state is exponentially small. Therefore, finding the")
    print("equilibrium values is equivalent to performing this entropy maximization.")

    print("\nSystem: We are considering a gas of photons, which are bosons. A key feature is that the number of photons")
    print("is not conserved. The state of the system can be described by the average number of photons, n_i, in each")
    print("available energy mode, ε_i.")
    
    # --- Define Symbolic Variables for the Equations ---
    n_i = sp.Symbol('n_i')        # Average occupation number of the i-th mode
    epsilon_i = sp.Symbol('ε_i')  # Energy of the i-th mode
    k_B = sp.Symbol('k_B')        # Boltzmann constant
    i = sp.Symbol('i')            # Index for summation
    E = sp.Symbol('E')            # Total mean energy
    S = sp.Symbol('S')            # Total entropy

    # --- 1. Equilibrium Mean Energy (E) ---
    print("\n--- Part 1: Equilibrium Mean Energy ---")
    print("The total mean energy of the system is the sum of the average energies in each mode.")
    print("The average energy in mode 'i' is the energy of that mode (ε_i) multiplied by the")
    print("average number of photons in it (n_i).")

    # Define the mathematical expression for total mean energy
    mean_energy_eq = sp.Eq(E, sp.Sum(n_i * epsilon_i, i))
    
    print("\nThe final equation for the equilibrium mean energy (E) is:")
    sp.pprint(mean_energy_eq, use_unicode=True)
    print("\nWhere:")
    print(f"  E:     Total mean energy")
    print(f"  n_i:   Average number of photons in the i-th energy state")
    print(f"  ε_i:   Energy of a single photon in the i-th state")
    print("  The sum is over all possible energy states 'i'.")

    # --- 2. Equilibrium Entropy (S) ---
    print("\n--- Part 2: Equilibrium Entropy ---")
    print("The entropy of a system of non-interacting bosons is derived from the number of possible")
    print("microscopic arrangements for a given macroscopic state {n_i}. Using combinatorial methods")
    print("or the partition function, we arrive at the following formula for the total entropy.")

    # Define the mathematical expression for the entropy term of a single mode
    # S_i/k_B = (1+n_i)ln(1+n_i) - n_i*ln(n_i)
    entropy_term = (1 + n_i) * sp.log(1 + n_i) - n_i * sp.log(n_i)
    
    # Define the equation for total entropy
    entropy_eq = sp.Eq(S, k_B * sp.Sum(entropy_term, i))
    
    print("\nThe final equation for the equilibrium entropy (S) is:")
    sp.pprint(entropy_eq, use_unicode=True)
    print("\nWhere:")
    print(f"  S:     Total entropy")
    print(f"  k_B:   Boltzmann constant")
    print(f"  n_i:   Average number of photons in the i-th energy state")
    print(f"  log:   Natural logarithm")
    print("  The sum is over all possible energy states 'i'.")

# Execute the function to print the solution
find_equilibrium_values()