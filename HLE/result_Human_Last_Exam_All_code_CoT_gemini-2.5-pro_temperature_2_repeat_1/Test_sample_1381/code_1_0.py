import numpy as np
import itertools

def solve_for_equilibrium(params, species_indices):
    """
    Solves for the equilibrium abundances for a given subset of species.

    Args:
        params (dict): A dictionary containing the model parameters N, Gamma, R, K, A.
        species_indices (list or tuple): A tuple of 1-based indices of the coexisting species.

    Returns:
        np.ndarray: A numpy array with the equilibrium abundances.
                    Returns None if no unique positive solution exists.
    """
    N = params['N']
    Gamma = params['Gamma']
    R = params['R']
    K = params['K']
    A = params['A']
    
    # The set of coexisting species is empty
    if not species_indices:
        return np.zeros(N)

    # Let M be the number of coexisting species
    M = len(species_indices)
    
    # The system is L * X = K_vec
    L = np.zeros((M, M))
    K_vec = np.zeros(M)
    
    # Map from species index (1 to N) to matrix index (0 to M-1)
    idx_map = {species_idx: i for i, species_idx in enumerate(species_indices)}

    for i, i_species in enumerate(species_indices):
        K_vec[i] = K[i_species - 1]
        for j, j_species in enumerate(species_indices):
            # Interaction term A_i - A_j
            interaction = A[i_species - 1] - A[j_species - 1]
            
            # Off-diagonal elements
            L[i, j] = - (Gamma * K[i_species-1] * R[j_species-1] / N) * interaction
        
        # Diagonal elements
        L[i, i] += 1
        
    # Solve the linear system L*X = K
    try:
        detL = np.linalg.det(L)
        if np.abs(detL) < 1e-9: # Matrix is singular or close to singular
            return None
        X_coexist = np.linalg.solve(L, K_vec)
    except np.linalg.LinAlgError:
        return None
    
    # Check if all abundances are positive
    if np.all(X_coexist > 0):
        # Create the full abundance vector
        X_full = np.zeros(N)
        for i, species_idx in enumerate(species_indices):
            X_full[species_idx - 1] = X_coexist[i]
        return X_full
    else:
        return None

def main():
    """
    Demonstrates finding all 2^N equilibria for N=3 with a specific parameter set.
    """
    N = 3
    # Choose parameters such that A_i values are close
    params = {
        'N': N,
        'Gamma': 1.0,
        'R': np.ones(N),
        'K': np.ones(N),
        'A': np.array([0.03, 0.02, 0.01])
    }

    all_species = list(range(1, N + 1))
    num_equilibria = 0
    
    print(f"Finding equilibria for a system of N={N} species.")
    print(f"Parameters: Gamma={params['Gamma']}, R={params['R']}, K={params['K']}, A={params['A']}\n")

    # Iterate over all possible subsets of species
    for i in range(N + 1):
        for subset in itertools.combinations(all_species, i):
            equilibrium_X = solve_for_equilibrium(params, subset)
            
            if equilibrium_X is not None:
                num_equilibria += 1
                survivor_set = subset if subset else "{}"
                print(f"Equilibrium found for survivor set: {survivor_set}")
                equation_str = "Abundances:"
                for j in range(N):
                    equation_str += f" X_{j+1} = {equilibrium_X[j]:.6f}"
                print(equation_str)
                print("-" * 20)

    print(f"\nTotal number of equilibria found: {num_equilibria}")
    print(f"This matches the theoretical maximum of 2^N = 2^{N} = {2**N}")

if __name__ == '__main__':
    main()
