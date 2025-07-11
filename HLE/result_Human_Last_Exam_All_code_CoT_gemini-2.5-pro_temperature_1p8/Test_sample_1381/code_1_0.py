import itertools

def analyze_equilibria(N):
    """
    This function analyzes and demonstrates the number of possible equilibria
    for the generalized Lotka-Volterra system with N species.

    The reasoning is as follows:
    The equilibrium conditions are given by a system of N quadratic equations.
    For each species i, the condition is X_i * P(X_1, ..., X_N) = 0, where P is a
    linear polynomial. This structure ensures that for each of the 2^N subsets of
    species, there is at most one equilibrium where precisely that subset has
    non-zero abundance. This gives an upper bound of 2^N equilibria.

    To show this maximum is achievable, we can select parameters such that all
    A_i values are equal. The system then simplifies to N uncoupled logistic
    growth equations. For each species i, its equilibrium abundance X_i* can be
    either 0 or its carrying capacity K_i.
    
    This function enumerates all 2^N equilibria for this simplified case.
    """
    
    print(f"Analysis for a system with N = {N} species.")
    print("---------------------------------------------")
    
    # In the simplified case (all A_i are equal), each species `i` can either
    # be extinct (X_i* = 0) or at its carrying capacity (X_i* = K_i).
    
    # Create symbolic representations for the carrying capacities.
    k_symbols = [f"K_{i+1}" for i in range(N)]
    
    # For each species, its state can be '0' or its corresponding 'K_i'.
    # We create a list of these state pairs for each species.
    species_states = [('0', k) for k in k_symbols]
    
    # Generate all possible equilibria by finding the Cartesian product of the states.
    # This corresponds to choosing one of the two states for each of the N species.
    equilibria_list = list(itertools.product(*species_states))
    
    num_equilibria = len(equilibria_list)
    
    print(f"By choosing parameters appropriately, we can obtain {num_equilibria} distinct equilibria.")
    print("These equilibria are of the form (X_1*, X_2*, ..., X_N*), where each X_i* is either 0 or K_i.")
    print("\nList of all possible equilibria points:")
    # Print the generated list of equilibria tuples
    for i, eq in enumerate(equilibria_list):
        print(f"{i+1}: {eq}")
        
    # The number of equilibria is calculated by the equation 2^N
    base = 2
    exponent = N
    final_equation_result = base**exponent
    
    print("\nThe number of equilibria is given by the final equation:")
    print(f"Number of equilibria = {base}**{exponent} = {final_equation_result}")

# Run the analysis for a specific number of species, N.
# The user can change this value to explore other cases.
N_species = 3
analyze_equilibria(N_species)