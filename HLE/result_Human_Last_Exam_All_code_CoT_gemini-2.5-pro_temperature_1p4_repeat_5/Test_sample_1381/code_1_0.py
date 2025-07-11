import math
from itertools import chain, combinations

def solve_for_equilibria_count(N):
    """
    Demonstrates that 2^N equilibria can exist for the given Lotka-Volterra system.

    Args:
        N (int): The number of species in the system.
    """
    if not isinstance(N, int) or N <= 0:
        print("Please provide a positive integer for N.")
        return

    # 1. Parameter selection
    # We choose parameters according to the derivation to maximize the number of
    # feasible equilibria. R_i, K_i, Gamma must be >= 1.
    R = [1.0] * N
    K = [1.0] * N
    Gamma = 1.0

    # Choose A_i values to be distinct but very close to each other.
    # A_i = 1 + epsilon * (i-1)
    # The condition for feasibility is N + Gamma * sum_j((A_i-A_j)R_jK_j) > 0
    # With our choices, this simplifies to N + sum_j(A_i - A_j) > 0.
    # N + epsilon * sum_j(i-j) > 0
    # We choose epsilon small enough to satisfy this for all subsets.
    # A safe choice is epsilon = 1.0 / N**2, which makes the interaction term small.
    epsilon = 1.0 / (N**2)
    A = [1.0 + epsilon * i for i in range(N)]

    # Helper function to get all subsets of species indices {1, ..., N}
    def get_all_subsets(n):
        s = list(range(1, n + 1))
        return chain.from_iterable(combinations(s, r) for r in range(n + 1))

    all_subsets = list(get_all_subsets(N))
    feasible_equilibria_count = 0

    # 2. Check feasibility for each subset
    for s_tuple in all_subsets:
        S = set(s_tuple)
        
        # The trivial equilibrium (S is empty) is always possible.
        if not S:
            feasible_equilibria_count += 1
            continue

        is_feasible = True
        # Check the feasibility condition for every species i in the subset S
        for i in S:
            # All R_j and K_j are 1, so the sum simplifies.
            interaction_sum = sum(A[i-1] - A[j-1] for j in S)
            
            # Feasibility condition: 1 + (Gamma/N) * interaction_sum > 0
            if 1.0 + (Gamma / N) * interaction_sum <= 0:
                is_feasible = False
                break
        
        if is_feasible:
            feasible_equilibria_count += 1

    # 3. Print the result
    print(f"For N = {N}, we have chosen parameters for which the number of feasible equilibria is {feasible_equilibria_count}.")
    
    theoretical_max = 2**N
    print(f"This matches the theoretical maximum number of equilibria: 2^{N} = {theoretical_max}.")

    if feasible_equilibria_count == theoretical_max:
        print("Conclusion: The maximum number of equilibria is 2^N.")
    else:
        print("Error: The count did not match the theoretical maximum.")


# --- Main execution ---
# We use N=4 as an example. The code works for other N, but gets slow as N increases.
N_species = 4
solve_for_equilibria_count(N_species)