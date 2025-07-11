import numpy as np
from scipy.optimize import fsolve
import warnings

# Suppress RuntimeWarning from log of a non-positive number
warnings.filterwarnings("ignore", category=RuntimeWarning)

def solve_integral_numerically():
    """
    This function provides a numerical example for the computation.
    It defines a set of system parameters, calculates the intermediate quantities μ₀ and Γ,
    and then computes the final resulting matrix's eigenvalues.
    """
    # --- Step 1: Define system parameters ---
    # We choose an undersampled case (M < N) where μ₀ < 0.
    M = 10
    N = 20.0 # Use float for division
    
    # Define the distinct eigenvalues (γ_r) and their multiplicities (K_r)
    # The sum of multiplicities must equal M.
    gamma_values = np.array([1.0, 2.5, 4.0])
    K_values = np.array([3, 5, 2])
    
    if np.sum(K_values) != M:
        raise ValueError("Sum of multiplicities K_r must be equal to M.")
        
    print("--- System Parameters ---")
    print(f"M = {M}, N = {N}")
    print(f"Eigenvalues of R (γ_r): {gamma_values}")
    print(f"Multiplicities (K_r): {K_values}")
    print("-------------------------\n")

    # --- Step 2: Numerically find μ₀ ---
    # μ₀ is a root of the equation: 1 - (1/N) * Σ [K_r * γ_r / (γ_r - μ)] = 0
    def mu0_equation(mu):
        return 1.0 - (1.0/N) * np.sum(K_values * gamma_values / (gamma_values - mu))

    # For the M < N case, the relevant root μ₀ is negative.
    # We provide an initial guess of -1.0 to the solver.
    mu0_solution = fsolve(mu0_equation, -1.0)
    mu0 = mu0_solution[0]
    
    print(f"--- Intermediate Calculations ---")
    print(f"Calculated μ₀ = {mu0:.8f}")

    # --- Step 3: Calculate Γ ---
    # Γ = (1/N) * Σ [K_r * (γ_r / (γ_r - μ₀))²]
    gamma_term = (gamma_values / (gamma_values - mu0))**2
    Gamma = (1.0/N) * np.sum(K_values * gamma_term)
    
    print(f"Calculated Γ = {Gamma:.8f}\n")

    # --- Step 4: Calculate the final expression ---
    # The value of the integral is the matrix -log(1 - Γ) * R * (R - μ₀I)⁻¹
    # This result is a diagonal matrix in the eigenbasis of R.
    # We compute its eigenvalues, which are given by:
    # λ'_r = -log(1 - Γ) * γ_r / (γ_r - μ₀)

    log_term = -np.log(1 - Gamma)
    result_eigenvalues = log_term * (gamma_values / (gamma_values - mu0))
    
    print("--- Final Result ---")
    print("The resulting matrix is diagonal in the same basis as R.")
    print("Its distinct eigenvalues are:")
    
    # "output each number in the final equation!"
    final_expression = f"Result Eigenvalues = -log(1 - {Gamma:.4f}) * [γ_r / (γ_r - ({mu0:.4f}))]"
    print(final_expression)
    
    for i, eig in enumerate(result_eigenvalues):
        print(f"For γ_{i+1} = {gamma_values[i]}, the resulting eigenvalue is: {eig:.8f}")
        
    # As the final answer requested by the user, we provide the simplified symbolic expression.
    final_answer_symbolic = "-log(1 - Γ) * R * (R - μ₀*I)⁻¹"
    
    # We output the components of this matrix, which are its eigenvalues.
    final_answer_numeric_list = [f"{eig:.8f}" for eig in result_eigenvalues]

if __name__ == '__main__':
    solve_integral_numerically()

# The symbolic expression for the value of the integral is:
# -log(1 - Γ) * R * (R - μ₀*I)⁻¹
# where Γ is given in the hint and μ₀ is the unique solution to 
# 1 = (1/N) * Σ[K_r * γ_r / (γ_r - μ₀)] such that the condition
# (1/N) * Σ[K_r * (γ_r / (γ_r - μ₀))²] < 1 is satisfied.
# The code above calculates the eigenvalues of this matrix for a specific numerical example.