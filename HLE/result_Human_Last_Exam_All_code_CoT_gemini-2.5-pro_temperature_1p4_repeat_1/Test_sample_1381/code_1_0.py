import numpy as np
import itertools
import sys

def count_equilibria(N):
    """
    Counts the number of feasible equilibria for a given N by constructing
    a suitable set of parameters.
    """
    # 1. Construct parameters R, K, Gamma, A such that all 2^N equilibria are feasible.
    # We choose R, K, Gamma to be simple values, and A values to be close to each other
    # to ensure the interaction terms are weak.
    R = np.ones(N)
    K = np.ones(N)
    Gamma = 1.0

    # We need A_i = Base + epsilon * i, where epsilon is small enough.
    # The condition is that 1 + (Gamma/N) * sum_{j in S} (A_i - A_j) > 0 for all S and i in S.
    # Let A_i = 2 + epsilon * i. Then A_i - A_j = epsilon * (i - j).
    # We need 1 + epsilon * (Gamma/N) * sum_{j in S} (i - j) > 0.
    # This requires epsilon < N / (Gamma * C_max) where C_max is the max of |sum(i-j)|
    C_max = 0
    indices = range(N)
    # Find the largest possible value of |sum(i-j)| to determine a safe epsilon.
    for i in indices:
        for s_size in range(1, N + 1):
            for s_indices in itertools.combinations(indices, s_size):
                if i in s_indices:
                    c_val = sum(i - j for j in s_indices)
                    if abs(c_val) > C_max:
                        C_max = abs(c_val)
    
    # Set epsilon to a safe, small value.
    if C_max > 0:
        epsilon = (N / Gamma) / (C_max * 1.1)
    else: # N=1 case
        epsilon = 0.1
    
    A = 2.0 + epsilon * np.arange(N)

    # 2. Iterate through all 2^N subsets S and check for feasible equilibria.
    equilibrium_count = 0
    
    for s_size in range(N + 1):
        for s_tuple in itertools.combinations(indices, s_size):
            S = list(s_tuple)
            
            # The all-zero equilibrium (for S = empty set) is always feasible.
            if not S:
                equilibrium_count += 1
                continue

            # For a non-empty subset S, solve for the equilibrium abundances.
            R_s, K_s, A_s = R[S], K[S], A[S]
            
            # Calculate coefficients for the linear system in U and V.
            S0 = np.sum(R_s * K_s)
            SA = np.sum(A_s * R_s * K_s)
            SA2 = np.sum(A_s**2 * R_s * K_s)
            
            # Solve the 2x2 linear system for U and V.
            G_over_N = Gamma / N
            m11 = 1 - G_over_N * SA
            m12 = G_over_N * S0
            m21 = -G_over_N * SA2
            m22 = 1 + G_over_N * SA
            
            M_uv = np.array([[m11, m12], [m21, m22]])
            C_uv = np.array([S0, SA])
            
            # This matrix is always invertible as per the analysis.
            U, V = np.linalg.solve(M_uv, C_uv)

            # Calculate the abundances X_i for i in S.
            X_s = K_s * (1 + G_over_N * A_s * U - G_over_N * V)
            
            # 3. A valid equilibrium requires all abundances to be positive.
            if np.all(X_s > 0):
                equilibrium_count += 1
    
    return equilibrium_count

if __name__ == '__main__':
    if len(sys.argv) > 1:
        try:
            N = int(sys.argv[1])
            if N < 1:
                raise ValueError
        except (ValueError, IndexError):
            print("Please provide a single positive integer N as a command-line argument.")
            print("Example: python your_script_name.py 3")
            sys.exit(1)
    else:
        # Default to N=4 if no argument is given.
        N = 4

    num_equilibria = count_equilibria(N)

    print(f"For N = {N}, a set of parameters can be chosen to allow for {num_equilibria} distinct equilibria.")
    print("This corresponds to the general formula 2^N, where N is the number of species.")
    print("\nThe final equation is:")
    # Using int() to ensure clean output without decimals
    print(f"{int(num_equilibria)} = 2^{N}")
