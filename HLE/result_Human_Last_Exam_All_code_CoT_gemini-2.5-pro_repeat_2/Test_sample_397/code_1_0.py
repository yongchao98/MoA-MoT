import numpy as np
from scipy.optimize import fsolve

def solve_integral_problem():
    """
    Solves for the value of the integral based on RMT results.

    The integral evaluates to -log(1 - Gamma) * I_M. This function
    computes the scalar part of this expression by first finding mu_0,
    then Gamma.
    """

    # --- Problem Parameters ---
    # These can be changed to analyze different systems.
    # M: dimension of observations
    # N: number of observations
    # gammas: distinct eigenvalues of the covariance matrix R
    # Ks: multiplicities of the corresponding eigenvalues
    M = 10
    N = 20
    gammas = np.array([2.0, 4.0])
    Ks = np.array([5, 5])
    
    # Ensure M is the sum of multiplicities
    if M != np.sum(Ks):
        raise ValueError("M must be the sum of all multiplicities Ks.")

    bar_M = len(gammas)
    print(f"Parameters: M={M}, N={N}, c=M/N={M/N:.2f}")
    print(f"Eigenvalues (gamma_r): {gammas}")
    print(f"Multiplicities (K_r):  {Ks}\n")


    # --- Step 1: Find mu_0 ---
    
    # Define the function whose roots (other than mu=0) we need.
    # 1 - (1/N) * sum(K_r * gamma_r / (gamma_r - mu)) = 0
    def mu_equation(mu):
        return 1.0 - (1/N) * np.sum(Ks * gammas / (gammas - mu))

    # Find the roots of the equation. We need to search for roots between the
    # poles of the function, which are the gammas.
    # The roots are located in the intervals (-inf, gamma_1), (gamma_1, gamma_2), ...
    
    # Sort gammas to define search intervals
    sorted_gammas = np.sort(gammas)
    
    potential_solutions = [0.0] # mu=0 is always a potential solution
    
    # Search for a root smaller than the smallest gamma
    # We search from a point far to the left of the smallest gamma
    search_start = sorted_gammas[0] - 100
    root, info, ier, msg = fsolve(mu_equation, x0=search_start, full_output=True)
    if ier == 1:
        potential_solutions.append(root[0])

    # Search for roots between consecutive gammas
    for i in range(bar_M - 1):
        # Starting point in the middle of the interval
        search_start = (sorted_gammas[i] + sorted_gammas[i+1]) / 2.0
        root, info, ier, msg = fsolve(mu_equation, x0=search_start, full_output=True)
        if ier == 1:
            potential_solutions.append(root[0])
            
    # Find the smallest solution
    mu_0 = np.min(potential_solutions)
    
    print("--- Intermediate Calculations ---")
    print(f"1. The smallest solution, mu_0, is: {mu_0:.8f}")

    # --- Step 2: Calculate Gamma ---
    
    gamma_ratio_sq = (gammas / (gammas - mu_0))**2
    Gamma = (1/N) * np.sum(Ks * gamma_ratio_sq)
    
    print(f"2. The value of Gamma is: {Gamma:.8f}")
    if Gamma >= 1:
        print("Warning: Gamma is >= 1, log(1-Gamma) is undefined for real numbers.")
        final_scalar_value = np.nan
    else:
        # --- Step 3: Compute the final value ---
        final_scalar_value = -np.log(1 - Gamma)

    print("\n--- Final Result ---")
    print("The integral evaluates to the matrix: -log(1 - Gamma) * I_M")
    print(f"The scalar value -log(1 - Gamma) is: {final_scalar_value:.8f}")
    print("\nThe final matrix is:")
    print(f"{final_scalar_value:.8f} * Identity({M})")


if __name__ == '__main__':
    solve_integral_problem()
