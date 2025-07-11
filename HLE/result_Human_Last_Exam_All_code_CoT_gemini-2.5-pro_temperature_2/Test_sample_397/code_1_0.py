import numpy as np
from scipy.optimize import fsolve

def solve_task():
    """
    Solves the complex integral problem by calculating its analytical solution.
    
    The steps are:
    1. Define the parameters M, N, and the eigenvalues/multiplicities of R.
    2. Numerically find the value of mu_0.
    3. Calculate Gamma using the provided formula.
    4. Construct the matrix R from its eigenvalues and projectors.
    5. Compute the final matrix Q(mu_0) = (R - mu_0 * I)^-1.
    6. Calculate the final result: (1 - Gamma) * Q(mu_0).
    7. Print the components of the final result equation.
    """
    # Step 1: Define parameters for a concrete example
    # Let's consider a simple case with M=3, N=5 (oversampled, c=0.6 < 1)
    M = 3
    N = 5
    
    # Distinct eigenvalues of R
    gamma_vals = np.array([2.0, 3.0, 5.0])
    # Multiplicities of eigenvalues
    K_vals = np.array([1, 1, 1])
    M_bar = len(gamma_vals)

    # Let's use the actual M for calculation, sum of K_vals should be M
    if np.sum(K_vals) != M:
        raise ValueError("Sum of multiplicities must be equal to M.")

    # Step 2: Find mu_0
    # mu_0 is a root of mu * (1 - (1/N) * sum(K_r * gamma_r / (gamma_r - mu))) = 0
    # The roots are mu=0 and the roots of the second term.
    # Let's define the function whose root we need
    def f_mu(mu):
        if np.any(np.isclose(mu, gamma_vals)):
            return float('inf')
        return 1.0 - (1.0/N) * np.sum(K_vals * gamma_vals / (gamma_vals - mu))

    # Based on the theory for c < 1, mu_0 = 0.
    # If M > N (c>1), there would be a negative root.
    if M < N:
        mu_0 = 0.0
    else: # M > N
        # For M > N, find the negative root of f_mu(mu) = 0
        # Search for a root in a negative interval
        mu_0_solution = fsolve(f_mu, -1.0)
        mu_0 = mu_0_solution[0]


    # Step 3: Calculate Gamma
    Gamma = (1.0/N) * np.sum(K_vals * (gamma_vals / (gamma_vals - mu_0))**2)

    # Step 4: Construct the matrix R
    # For simplicity, let R be a diagonal matrix with the given eigenvalues.
    # The projectors P_m are diagonal matrices with 1 at the corresponding position.
    R = np.diag(gamma_vals)
    
    # Step 5: Calculate Q(mu_0)
    I_M = np.identity(M)
    # The term (R - mu_0*I) is invertible since mu_0 is not an eigenvalue of R.
    Q_mu0 = np.linalg.inv(R - mu_0 * I_M)
    
    # Step 6: Calculate the final result
    result_matrix = (1 - Gamma) * Q_mu0

    # Step 7: Print the components of the equation and the final result
    print("Problem Parameters:")
    print(f"M = {M}, N = {N}")
    print(f"Eigenvalues of R (gamma_r): {gamma_vals}")
    print(f"Multiplicities (K_r): {K_vals}\n")
    
    print("Intermediate Calculations:")
    print(f"Calculated mu_0 = {mu_0:.4f}")
    print(f"Calculated Gamma = {Gamma:.4f}")
    print(f"1 - Gamma = {1 - Gamma:.4f}\n")

    print("Matrices:")
    print("R =")
    print(R)
    print("\nQ(mu_0) = (R - mu_0*I)^-1 =")
    print(Q_mu0)
    
    print("\nFinal Result:")
    print("The integral evaluates to (1 - Gamma) * Q(mu_0):")
    print(f"({1-Gamma:.4f}) * Q(mu_0) =")
    print(result_matrix)
    print(f"\n<<< {result_matrix} >>>")


solve_task()
