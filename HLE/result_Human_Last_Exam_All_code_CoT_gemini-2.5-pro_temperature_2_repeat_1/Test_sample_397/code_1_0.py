import numpy as np
from scipy.optimize import fsolve
import math

def calculate_integral_value():
    """
    This function implements the derived analytical solution for the integral.
    It sets up a numerical example, calculates all the intermediate quantities
    (mu_0, Gamma), and computes the final matrix value of the integral.
    """
    # Plan:
    # 1. Set up the parameters for the problem: M, N, and the eigenvalues of R.
    #    We will use an example in the undersampled case (M > N), where the result is non-trivial.
    # 2. Define the true covariance matrix R based on the given eigenvalues.
    # 3. Find the value of mu_0 by numerically solving the equation that defines it.
    #    mu_0 is the unique negative solution to N - M = mu * sum(K_r / (gamma_r - mu)).
    # 4. Calculate Gamma using the derived value of mu_0.
    # 5. Compute the final matrix expression: I = -mu_0 * inv(R - mu_0*I) * log(1-Gamma).
    # 6. Print all the components of the final formula and the final result.

    # Step 1: Define parameters
    M = 4  # Dimension of observations
    N = 2  # Number of observations
    # M > N corresponds to the undersampled case where mu_0 < 0.

    # Eigenvalues of R and their multiplicities
    gammas = np.array([1.0, 3.0])
    Ks = np.array([2, 2])
    if np.sum(Ks) != M:
        raise ValueError("Sum of multiplicities must be equal to M.")

    # Step 2: Construct the matrix R
    # R is represented by a diagonal matrix with its eigenvalues.
    R_diag = []
    for gamma, K in zip(gammas, Ks):
        R_diag.extend([gamma] * K)
    R = np.diag(R_diag)
    I_M = np.identity(M)

    # Step 3: Find mu_0
    # Define the function whose root defines mu_0 for the undersampled case.
    # f(mu) = mu * sum(K_r / (gamma_r - mu)) - (N - M) = 0
    def equation_for_mu(mu):
        if mu >= np.min(gammas): # Avoid poles at eigenvalues
            return float('inf')
        sum_term = np.sum([K / (gamma - mu) for gamma, K in zip(gammas, Ks)])
        return mu * sum_term - (N - M)

    # Find the negative root of the equation using a numerical solver.
    # For mu < 0, the function is monotonic, so a root is guaranteed.
    mu_0_solution = fsolve(equation_for_mu, x0=-0.1)
    mu_0 = mu_0_solution[0]

    # Step 4: Calculate Gamma
    # Gamma = (1/N) * sum(K_r * (gamma_r / (gamma_r - mu_0))^2)
    gamma_sum_term = np.sum([K * (gamma / (gamma - mu_0))**2 for gamma, K in zip(gammas, Ks)])
    Gamma = (1/N) * gamma_sum_term

    if not 1 - Gamma > 0:
        raise ValueError("Condition 1 - Gamma > 0 is not met.")

    # Step 5: Compute the final integral value based on the analytical formula
    log_term = np.log(1 - Gamma)
    matrix_term = -mu_0 * np.linalg.inv(R - mu_0 * I_M)
    integral_value = matrix_term * log_term

    # Step 6: Print the components and the final result as requested
    print("--- Calculating the value of the integral ---")
    print("The analytical formula is: I = -mu_0 * inv(R - mu_0*I) * log(1-Gamma)")
    print("\n--- Component Values ---")
    print(f"M = {M}")
    print(f"N = {N}")
    print(f"\nEigenvalues of R (gamma_r): {gammas}")
    print(f"Multiplicities (K_r): {Ks}")
    print("\nThe matrix R:")
    print(R)
    print(f"\nCalculated mu_0 = {mu_0:.4f}")
    print(f"Calculated Gamma = {Gamma:.4f}")

    print("\n--- Final Equation with Numbers ---")
    print("Scalar part: log(1 - Gamma) = log(1 - {:.4f}) = {:.4f}".format(Gamma, log_term))
    print("\nMatrix part: -mu_0 * inv(R - mu_0*I) = ")
    print(matrix_term)

    print("\n--- Final Value of the Integral (I) ---")
    print(integral_value)

calculate_integral_value()