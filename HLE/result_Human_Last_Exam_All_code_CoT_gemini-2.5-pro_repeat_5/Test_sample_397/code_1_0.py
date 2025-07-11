import numpy as np
from scipy.optimize import brentq

def solve_integral_value():
    """
    Calculates the scalar value of the complex integral based on the provided RMT definitions.
    
    The final value is derived to be (1 - M/N) * log(1 - Gamma) for M < N and 0 for M >= N.
    """
    # Problem parameters
    # M-dimensional observations
    M = 50
    # N observations
    N = 100
    # Distinct eigenvalues of R
    gammas = np.array([1.0, 2.5, 4.0])
    # Multiplicities of the eigenvalues
    Ks = np.array([20, 20, 10])
    
    # Check for consistency
    if np.sum(Ks) != M:
        raise ValueError("Sum of multiplicities Ks must be equal to M.")

    print(f"M = {M}, N = {N}")
    print("Eigenvalues (gamma_r):", gammas)
    print("Multiplicities (K_r):", Ks)
    print("-" * 30)

    c = M / N
    print(f"Ratio c = M/N = {c:.4f}")

    if c >= 1:
        # Oversampled or critically sampled case
        # In this regime, mu_0 = 0, which makes the resulting integral value 0.
        print("Oversampled or critical regime (c >= 1).")
        result = 0.0
        mu_0 = 0.0
        gamma_val = M / N # For c > 1, mu_0=0, Gamma = M/N
        print(f"mu_0 = {mu_0}")
        print(f"Gamma = {gamma_val:.4f}")
        print(f"Final scalar value of the integral is {result}")
        print("\nFinal equation:")
        print(f"Result = (1 - {c:.4f}) * log(1 - {gamma_val:.4f}) = {result}")

    else:
        # Undersampled case
        print("Undersampled regime (c < 1).")
        
        # Define the function to find the root mu_0
        # The equation is: 1 - c = (mu/N) * sum(K_r / (gamma_r - mu))
        def mu_equation(mu):
            if mu in gammas:
                return np.inf
            return 1 - c - (mu / N) * np.sum(Ks / (gammas - mu))

        # The solution mu_0 is in (0, min(gammas))
        gamma_1 = np.min(gammas)
        try:
            mu_0 = brentq(mu_equation, 0, gamma_1 * 0.9999)
        except ValueError:
            print("Could not find a root for mu_0 in the interval (0, gamma_1).")
            return

        print(f"Numerically found mu_0 = {mu_0:.4f}")

        # Calculate Gamma using mu_0
        gamma_val = (1 / N) * np.sum(Ks * (gammas / (gammas - mu_0))**2)
        print(f"Gamma = {gamma_val:.4f}")

        # The final result is (1 - c) * log(1 - Gamma)
        # We must have Gamma < 1 for the log to be real.
        if gamma_val >= 1:
            print("Warning: Gamma is >= 1, log will be undefined or complex.")
            # For z' > 0 condition, we must have 1-Gamma > 0
            result = np.nan
        else:
            result = (1 - c) * np.log(1 - gamma_val)

        print(f"Final scalar value of the integral is {result:.4f}")
        
        # Output the final equation with all the numbers
        print("\nFinal equation:")
        print(f"Result = (1 - {M}/{N}) * log(1 - {gamma_val:.4f})")
        print(f"Result = ({1-c:.4f}) * log({1-gamma_val:.4f})")
        print(f"Result = ({1-c:.4f}) * ({np.log(1 - gamma_val):.4f}) = {result:.4f}")

solve_integral_value()