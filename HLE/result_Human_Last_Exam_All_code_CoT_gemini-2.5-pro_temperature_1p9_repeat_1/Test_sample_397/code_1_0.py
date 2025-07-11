import numpy as np
import sympy as sp

def solve_integral():
    """
    Solves the integral for a specific case and prints the result.
    """
    # Define parameters of the problem
    M = 2
    N = 4
    gammas = np.array([3.0, 5.0])
    Ks = np.array([1, 1])
    R_diag = np.array([3.0, 5.0])
    R = np.diag(R_diag)
    
    # Define the symbolic variable for mu
    mu = sp.Symbol('mu')
    
    # Set up the equation for mu_k
    # mu * (1 - (1/N) * sum(K_r * gamma_r / (gamma_r - mu))) = 0
    sum_term = sum(K * g / (g - mu) for g, K in zip(gammas, Ks))
    equation = mu * (1 - (1/N) * sum_term)
    
    # Solve for mu_k. We look for real solutions.
    solutions = [s for s in sp.solve(equation, mu) if s.is_real]
    
    # Identify mu_0. Based on convention, it is usually the smallest root.
    # The problem statement has some ambiguity, so we choose one for demonstration.
    # Let's pick the smallest non-zero solution for this example.
    # In a fully specified problem, mu_0 would be uniquely identified.
    mu_k = sorted([float(s) for s in solutions])
    if len(mu_k) < len(gammas) + 1:
        raise ValueError("Could not find all expected roots for mu.")
        
    mu_0 = mu_k[1] # Choose mu_0=1 based on sorted solutions [0, 1, 7]

    # Calculate Gamma
    gamma_val = (1/N) * np.sum(Ks * (gammas / (gammas - mu_0))**2)
    
    # Calculate the matrix R - mu_0*I and its inverse
    I = np.identity(M)
    R_minus_mu0I = R - mu_0 * I
    Q_mu0 = np.linalg.inv(R_minus_mu0I)

    # Calculate the logarithm term, which might be complex
    log_term = np.log(1 - gamma_val)

    # Calculate the final result matrix
    result_matrix = -mu_0 * log_term * Q_mu0

    # Output the results
    print("Assumed Parameters:")
    print(f"M = {M}, N = {N}")
    print(f"Eigenvalues (gamma_r): {gammas}")
    print(f"Multiplicities (K_r): {Ks}")
    print("-" * 30)
    print("Intermediate Calculations:")
    print(f"Solutions for mu are: {[round(s, 4) for s in mu_k]}")
    print(f"Chosen mu_0 = {mu_0:.4f}")
    print(f"Gamma = {gamma_val:.4f}")
    print(f"log(1 - Gamma) = {log_term:.4f}")
    print("-" * 30)
    print("Final Result:")
    print("The value of the integral is a matrix given by the expression:")
    print("-mu_0 * log(1 - Gamma) * (R - mu_0*I)^(-1)")
    print("\nCalculated value for mu_0:")
    print(f"mu_0 = {mu_0:.4f}")
    print("\nCalculated value for log(1 - Gamma):")
    print(f"log(1 - Gamma) = {log_term.real:.4f} + {log_term.imag:.4f}j")
    print("\nCalculated matrix (R - mu_0*I)^(-1):")
    print(Q_mu0)
    print("\nFinal Result Matrix:")
    print(result_matrix)
    print("\nIn equation form, for each diagonal element i:")
    final_diag = np.diag(result_matrix)
    print(f"Result_ii = -({mu_0:.4f}) * ({log_term:.4f}) * (1 / ({R_diag[0]:.4f} - {mu_0:.4f})) = {final_diag[0]:.4f}")
    print(f"Result_ii = -({mu_0:.4f}) * ({log_term:.4f}) * (1 / ({R_diag[1]:.4f} - {mu_0:.4f})) = {final_diag[1]:.4f}")

solve_integral()