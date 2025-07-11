import numpy as np
from scipy.optimize import fsolve

def solve_integral_numerically(M, N, gammas, Ks):
    """
    Calculates the value of the complex integral based on the derived non-zero result.

    Args:
        M (int): Dimension of the observation space.
        N (int): Number of observations.
        gammas (list or np.array): List of distinct eigenvalues of R.
        Ks (list or np.array): List of multiplicities for each eigenvalue.
    """
    if M != sum(Ks):
        raise ValueError("Sum of multiplicities Ks must be equal to M.")
    if len(gammas) != len(Ks):
        raise ValueError("Length of gammas and Ks must be the same.")

    # Construct the full R matrix for the final calculation
    # In the eigenbasis, R is diagonal.
    R_diag = np.repeat(gammas, Ks)
    R = np.diag(R_diag)

    # Define the equation for mu_0
    # 1 - (1/N) * sum(K_r * gamma_r / (gamma_r - mu)) = 0
    def mu_equation(mu):
        return 1.0 - (1.0 / N) * np.sum(Ks * gammas / (gammas - mu))

    # Find mu_0. It's the root outside the interval [min(gammas), max(gammas)].
    # We search for a root starting from a point smaller than the smallest eigenvalue.
    mu_0_initial_guess = np.min(gammas) - 1 
    mu_0 = fsolve(mu_equation, mu_0_initial_guess)[0]

    # Calculate Gamma
    gamma_term = gammas / (gammas - mu_0)
    Gamma = (1.0 / N) * np.sum(Ks * (gamma_term**2))

    # Calculate the matrix Q(mu_0) = (R - mu_0 * I)^-1
    Q_mu_0 = np.linalg.inv(R - mu_0 * np.identity(M))

    # Calculate the final expression for the integral
    # I = -mu_0 * log(1 - Gamma) * Q(mu_0)
    # Note: 1 - Gamma can be negative, so we use np.log for complex numbers.
    log_term = np.log(1 - Gamma)
    integral_value = -mu_0 * log_term * Q_mu_0
    
    # Print the final equation with computed values
    print("Based on the interpretation that the contour encloses the pole, the integral evaluates to:")
    print(f"I = -mu_0 * log(1 - Gamma) * (R - mu_0*I)^-1")
    print("\nWhere the computed values are:")
    print(f"mu_0 = {mu_0}")
    print(f"Gamma = {Gamma}")
    print(f"log(1 - Gamma) = {log_term}")
    
    print("\nThe final matrix result is:")
    # To satisfy the prompt "output each number in the final equation!",
    # we print the components of the final result.
    # The equation is I = C * A, where C is a scalar and A is a matrix.
    scalar_C = -mu_0 * log_term
    matrix_A = Q_mu_0
    print(f"I = ({scalar_C}) * A, where A is the matrix:")
    print(matrix_A)
    print("\nResulting in I =")
    print(integral_value)


# Example usage with some parameters
M_val = 2
N_val = 4
# Distinct eigenvalues of R
gammas_val = np.array([1.0, 2.0])
# Multiplicities of eigenvalues
Ks_val = np.array([1, 1])

solve_integral_numerically(M_val, N_val, gammas_val, Ks_val)

<<<
I = -mu_0 * log(1 - Gamma) * (R - mu_0*I)^-1
>>>