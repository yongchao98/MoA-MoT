import itertools

def find_maximum_equilibria(N, K_values):
    """
    This function demonstrates the case for the maximum number of equilibria
    in the generalized Lotka-Volterra system.

    The maximum number of equilibria is 2^N. This is achieved when all
    parameters A_i are equal, which simplifies the system to N independent
    logistic equations. In this case, the equilibria are all possible
    combinations where each species i is either at abundance 0 or at its
    carrying capacity K_i.

    Args:
        N (int): The number of species.
        K_values (list): A list of carrying capacities [K_1, K_2, ..., K_N].

    Returns:
        list: A list of all possible equilibrium points (tuples).
    """
    if not isinstance(N, int) or N < 1:
        raise ValueError("N must be a positive integer.")
    if len(K_values) != N:
        raise ValueError(f"The length of K_values ({len(K_values)}) must be equal to N ({N}).")

    # For each species i, its equilibrium abundance can be 0 or K_i.
    # We create a list of these possible pairs for all N species.
    possible_abundances_per_species = []
    for k in K_values:
        possible_abundances_per_species.append((0, k))

    # itertools.product generates the Cartesian product of the input iterables.
    # This gives all possible combinations of abundances.
    equilibria = list(itertools.product(*possible_abundances_per_species))
    
    return equilibria

def main():
    """
    Main function to execute the analysis for a given N.
    """
    # Let's consider a system with N species.
    # You can change this value to see the result for different numbers of species.
    N = 3

    # We need to define the carrying capacities K_i for each species.
    # The parameters R_i and Gamma must be >= 1.
    # The crucial choice is to set all A_i to be equal (and >= 1).
    # For this example, we assume R_i=1, Gamma=1, A_i=1 for all i.
    # Let's set some example carrying capacities.
    K = [100, 150, 200]
    
    # Calculate the equilibria for this specific case
    all_equilibria = find_maximum_equilibria(N, K)
    num_equilibria = len(all_equilibria)

    print(f"For a system of N = {N} species:")
    print("The maximum number of equilibria is given by the equation: 2^N.")
    # The instruction "output each number in the final equation" is interpreted here.
    print(f"Calculation: 2^{N} = {num_equilibria}")
    print("\nThis maximum is achieved when parameters are chosen such that the interaction terms cancel out")
    print("(e.g., by setting all A_i to be equal).")
    
    print(f"\nFor N = {N} and carrying capacities K = {K}, the {num_equilibria} possible equilibria are:")
    for i, eq in enumerate(all_equilibria):
        print(f"{i+1:2d}: {eq}")

if __name__ == "__main__":
    main()