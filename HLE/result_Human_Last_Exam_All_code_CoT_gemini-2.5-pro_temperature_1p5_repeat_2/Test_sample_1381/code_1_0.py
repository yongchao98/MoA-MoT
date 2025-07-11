import numpy as np
import itertools

def find_max_equilibria(N):
    """
    Demonstrates that 2**N equilibria can exist for the given Lotka-Volterra system.

    Args:
        N (int): The number of species in the system.
    """
    if N > 15:
        print(f"N={N} is too large to check all 2**N subsets in a reasonable time.")
        return

    # 1. Define the system parameters.
    # We choose R, K, A, and a very small Gamma to approximate the uncoupled case.
    R = np.full(N, 2.0)
    K = np.full(N, 10.0)
    A = np.arange(1.0, N + 1.0)
    Gamma = 1e-6

    # Counter for valid equilibria found
    equilibrium_count = 0
    
    # 2. Iterate over all possible subsets of species that could be present.
    # An equilibrium is defined by the subset of species with non-zero abundance.
    # The total number of subsets of {0, ..., N-1} is 2**N.
    indices = range(N)
    for s_size in range(N + 1):
        for s_indices in itertools.combinations(indices, s_size):
            # s_indices is a tuple representing the subset S of present species.
            S = list(s_indices)
            
            # The case S = {} corresponds to the trivial equilibrium where all species are extinct.
            if not S:
                equilibrium_count += 1
                continue

            # 3. For each non-empty subset S, solve for the equilibrium abundances.
            s = len(S)
            M_S = np.zeros((s, s))
            R_S = np.array([R[i] for i in S])

            for i_idx, i in enumerate(S):
                for j_idx, j in enumerate(S):
                    if i == j:
                        # Diagonal elements of the interaction matrix
                        M_S[i_idx, j_idx] = R[i] / K[i]
                    else:
                        # Off-diagonal elements
                        interaction = (Gamma / N) * R[i] * R[j] * (A[i] - A[j])
                        M_S[i_idx, j_idx] = -interaction
            
            # 4. Solve the linear system M_S * X_S = R_S
            try:
                # np.linalg.solve finds X_S for the system
                X_S = np.linalg.solve(M_S, R_S)
                
                # An equilibrium is physically valid only if all present species have positive abundance.
                if np.all(X_S > 1e-9): # Use a small tolerance for floating point comparison
                    equilibrium_count += 1
            except np.linalg.LinAlgError:
                # If the matrix M_S is singular, no unique solution exists for this subset.
                # For generic parameters (especially small Gamma), this is unlikely.
                pass

    # 5. Print the results and the final equation.
    print(f"For a system with N = {N} species:")
    print(f"Number of valid equilibria found by numerical simulation: {equilibrium_count}")
    theoretical_max = 2**N
    print(f"Theoretical maximum number of equilibria is 2^N.")
    print("\nFinal equation for the maximum number of equilibria:")
    print(f"count = 2^{N} = {theoretical_max}")


if __name__ == '__main__':
    # You can change N here to test for different numbers of species.
    # The calculation time grows as 2^N. N=4 is quick.
    N_species = 4
    find_max_equilibria(N_species)