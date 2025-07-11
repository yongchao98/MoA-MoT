import numpy as np
from scipy.linalg import eigh

def solve_eigenvalues(n, k, N=2000, R=30):
    """
    Solves for the eigenvalues of the radial operator L_k for a given n and k.
    
    n: dimension of the catenoid (surface is in R^{n+1})
    k: spherical harmonic index
    N: number of grid points
    R: domain limit [-R, R]
    """
    # Discretize the domain, avoiding the endpoints for Dirichlet boundary conditions
    rho = np.linspace(-R, R, N + 2)
    rho_internal = rho[1:-1]
    d_rho = rho[1] - rho[0]

    # Define the functions from the operator expression
    # Use a small epsilon to avoid division by zero at rho=0
    eps = 1e-9
    
    def F_rho_abs(r):
        # Regularized version of |F_rho| to handle the singularity at rho=0
        # The numerator has abs(r)
        # The term under the square root is approximated near r=0 to avoid nan.
        num = np.abs(r) * (r**2 + 1)**((n - 2) / 2)
        den_sq = (r**2 + 1)**(n - 1) - 1
        # Use Taylor expansion for small r to ensure stability
        den = np.sqrt(np.maximum(den_sq, (n - 1) * r**2))
        return num / (den + eps)

    # w and p functions for the Sturm-Liouville operator
    w = (rho_internal**2 + 1)**((n - 1) / 2) * F_rho_abs(rho_internal)
    p = (rho_internal**2 + 1)**((n - 1) / 2) / (F_rho_abs(rho_internal) + eps)
    
    # Potential function V_k
    lambda_k = k * (k + n - 2)
    V_k = -lambda_k / (rho_internal**2 + 1) + n * (n - 1) / (rho_internal**2 + 1)**n

    # Assemble the matrices for the generalized eigenvalue problem A*f = lambda*B*f
    # B is the diagonal weight matrix
    B = np.diag(w)
    
    # A is the symmetric matrix from the differential operator and potential
    p_half = (p[:-1] + p[1:]) / 2
    
    main_diag = - (p_half[:-1] + p_half[1:]) / d_rho**2 + V_k[1:-1] * w[1:-1]
    off_diag = p_half[1:-1] / d_rho**2
    
    A = np.diag(main_diag) + np.diag(off_diag, k=1) + np.diag(off_diag, k=-1)

    # Solve the generalized eigenvalue problem
    # We only need the eigenvalues
    eigenvalues = eigh(A, B[1:-1, 1:-1], eigvals_only=True)
    
    # Count positive eigenvalues (with a small tolerance)
    positive_eigenvalue_count = np.sum(eigenvalues > 1e-6)
    
    return positive_eigenvalue_count

def main():
    """
    Main function to calculate the number of positive eigenvalues.
    """
    print("Finding the number of positive eigenvalues for the stability operator L.")
    print("The operator is given by:")
    print("L = (1/(<ρ>^(n-1)|F_ρ|)) * d/dρ(<ρ>^(n-1)|F_ρ|⁻¹ * d/dρ) + (1/<ρ>²) * Δ_S + n(n-1)/<ρ>^(2n)")
    print("where <ρ> = sqrt(ρ² + 1) and |F_ρ| is a function of ρ and n.")
    print("\nWe will solve this numerically for a representative dimension n=3.")
    
    n = 3
    print(f"Set n = {n}")
    
    total_positive_eigenvalues = 0
    
    # We check for k = 0, 1, 2. For k>2, positive eigenvalues are not expected.
    for k in range(3):
        count = solve_eigenvalues(n=n, k=k)
        print(f"For spherical harmonic index k = {k}, found {count} positive eigenvalue(s).")
        total_positive_eigenvalues += count
        
    print("\n" + "="*40)
    print(f"The total number of positive eigenvalues for the operator L is: {total_positive_eigenvalues}")
    print("="*40)

if __name__ == "__main__":
    main()