import numpy as np
from scipy.optimize import brentq

def solve_integral_numerically():
    """
    This script calculates the value of the specified complex integral for a
    given set of parameters.

    The analytical result of the integral is:
    J = -mu_0 * (R - mu_0*I)^-1 * log(1 - Gamma)

    where:
    - R is the population covariance matrix.
    - mu_0 is the smallest real root of z(mu) = 0.
    - Gamma = (1/N) * sum(K_r * (gamma_r / (gamma_r - mu_0))^2)

    The script defines a sample problem, calculates mu_0 numerically,
    computes the necessary components, and prints the final result matrix.
    """

    # --- Problem Parameters ---
    # We choose an undersampled case (M > N) for a non-zero result.
    M = 10  # Dimension of observations
    N = 5   # Number of observations
    
    # Distinct eigenvalues of R and their multiplicities
    gamma_r = np.array([1.0, 2.0, 3.0])
    K_r = np.array([3, 4, 3])
    
    # Sanity check for dimensions
    if np.sum(K_r) != M:
        raise ValueError("Sum of multiplicities K_r must be equal to M.")

    print("--- Problem Setup ---")
    print(f"M = {M}, N = {N}")
    print(f"Eigenvalues gamma_r = {gamma_r}")
    print(f"Multiplicities K_r = {K_r}\n")

    # --- Step 1: Find mu_0 ---
    # mu_0 is the smallest real root of z(mu)=0, which for mu != 0 is equivalent to
    # f(mu) = 1 - (1/N) * sum(K_r * gamma_r / (gamma_r - mu)) = 0.
    # For the undersampled case M > N, we expect mu_0 < 0.
    def find_mu0_func(mu):
        return 1.0 - (1/N) * np.sum(K_r * gamma_r / (gamma_r - mu))

    # We search for the root in a large negative interval.
    # The upper bound is the smallest eigenvalue, min(gamma_r).
    try:
        mu_0 = brentq(find_mu0_func, -100, np.min(gamma_r) - 1e-9)
    except ValueError:
        print("Could not find a root for mu_0 in the search interval.")
        return

    print("--- Intermediate Calculations ---")
    print(f"Found mu_0 = {mu_0:.8f}")

    # --- Step 2: Construct Matrix R and calculate Gamma ---
    # R is a diagonal matrix with eigenvalues gamma_r
    eigenvalues = np.repeat(gamma_r, K_r)
    R = np.diag(eigenvalues)

    # Calculate Gamma
    gamma_term = (gamma_r / (gamma_r - mu_0))**2
    Gamma = (1/N) * np.sum(K_r * gamma_term)
    print(f"Calculated Gamma = {Gamma:.8f}")

    log_term = np.log(1 - Gamma)
    print(f"Calculated log(1 - Gamma) = {log_term:.8f}\n")
    if 1 - Gamma <= 0:
        print("Warning: log argument (1-Gamma) is not positive. The result might be complex or undefined.")
    
    # --- Step 3: Calculate the final result matrix J ---
    # J = -mu_0 * log(1-Gamma) * (R - mu_0*I)^-1
    # Since R is diagonal, (R - mu_0*I)^-1 is also diagonal.
    Q_mu0_diag = 1.0 / (eigenvalues - mu_0)
    Q_mu0 = np.diag(Q_mu0_diag)
    
    J = -mu_0 * log_term * Q_mu0

    # --- Step 4: Output the results ---
    print("--- Final Result: The value of the integral is the matrix J ---")
    # Since J is diagonal, we print its diagonal for brevity
    print("The resulting matrix J is diagonal. Its diagonal elements are:")
    final_equation_template = "-({:.4f}) * ({:.4f}) * (1 / ({} - ({:.4f})))"

    for i in range(M):
      gamma_i = eigenvalues[i]
      j_ii = J[i, i]
      final_equation = final_equation_template.format(mu_0, log_term, gamma_i, mu_0)
      print(f"J[{i},{i}] = {final_equation} = {j_ii:.8f}")


# Execute the function
solve_integral_numerically()