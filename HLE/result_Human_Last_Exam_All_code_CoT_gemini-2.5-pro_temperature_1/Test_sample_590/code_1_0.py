import numpy as np
from scipy.linalg import eig

def solve_and_count():
    """
    Numerically solves for the number of positive eigenvalues of the catenoid
    stability operator for the n=2 case.
    """
    n = 2  # Dimension of the catenoid
    N = 1000  # Number of grid points
    R = 20.0  # Domain for rho is [-R, R]

    rho, h = np.linspace(-R, R, N, retstep=True)

    # Finite difference matrix for the second derivative (Dirichlet boundary conditions)
    D2 = (np.diag(np.ones(N - 1), 1) - 2 * np.diag(np.ones(N)) + np.diag(np.ones(N - 1), -1)) / h**2

    # Finite difference matrix for the first derivative
    D1 = (np.diag(np.ones(N - 1), 1) - np.diag(np.ones(N - 1), -1)) / (2 * h)
    # Use forward/backward difference at the boundaries
    D1[0, 0] = -1.5/h; D1[0, 1] = 2/h; D1[0, 2] = -0.5/h
    D1[N-1, N-1] = 1.5/h; D1[N-1, N-2] = -2/h; D1[N-1, N-3] = 0.5/h
    
    ab_rho = np.sqrt(rho**2 + 1)
    
    # The radial part of the Laplace-Beltrami operator on the catenoid
    # For n=2, this is d^2/drho^2 + (rho/(rho^2+1)) * d/drho
    laplacian_rho_matrix = D2 + np.diag(rho / ab_rho**2) @ D1

    # --- Count positive eigenvalues for k=0 mode ---
    k0 = 0
    # Potential V_0 = 2/cosh^4(rho) = 2/angle_bracket(rho)^4
    potential_V0 = 2.0 / ab_rho**4
    L0_matrix = laplacian_rho_matrix + np.diag(potential_V0)
    eigenvalues_k0 = np.real(eig(L0_matrix, right=False))
    # Count eigenvalues > 0, with a small tolerance for numerical precision
    count_k0 = np.sum(eigenvalues_k0 > 1e-6)

    # --- Count positive eigenvalues for k=1 mode ---
    k1 = 1
    # Eigenvalue of Laplacian on S^1 is -k^2 = -1
    # Potential V_1 = 2/cosh^4(rho) - 1/cosh^2(rho)
    potential_V1 = (2.0 / ab_rho**4) - (k1**2 / ab_rho**2)
    L1_matrix = laplacian_rho_matrix + np.diag(potential_V1)
    eigenvalues_k1 = np.real(eig(L1_matrix, right=False))
    count_k1 = np.sum(eigenvalues_k1 > 1e-6)

    # Total number of positive eigenvalues = sum of counts * multiplicity
    # Multiplicity of k=0 on S^1 is 1
    # Multiplicity of k=1 on S^1 is 2
    multiplicity_k0 = 1
    multiplicity_k1 = 2
    total_pos_eigenvalues = multiplicity_k0 * count_k0 + multiplicity_k1 * count_k1

    # The final output should be the equation showing the calculation
    print(f"{multiplicity_k0} * {count_k0} + {multiplicity_k1} * {count_k1} = {total_pos_eigenvalues}")


solve_and_count()