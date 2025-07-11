import numpy as np
from scipy.linalg import eig

def count_positive_eigenvalues(n, k, num_points=2000, rho_max=20):
    """
    Counts the number of positive eigenvalues for the given operator L_k for a given n.

    Args:
        n (int): Dimension parameter from the problem.
        k (int): Spherical harmonic mode number.
        num_points (int): Number of points for discretization.
        rho_max (float): Maximum value of rho for domain truncation.

    Returns:
        int: The number of positive eigenvalues found.
    """
    rho = np.linspace(-rho_max, rho_max, num_points)
    h = rho[1] - rho[0]

    # Define helper functions
    def japanese_bracket(r):
        return np.sqrt(r**2 + 1)

    def F_rho_abs(r, n_val):
        # Handle rho=0 case to avoid division by zero
        r_safe = np.where(r == 0, 1e-9, r)
        jb_val = japanese_bracket(r_safe)
        
        # Use a stable calculation for the denominator
        # For small rho, (rho^2+1)^(n-1) - 1 is approx (n-1)rho^2
        if n_val == 2:
            # For n=2, |F_rho| simplifies to 1 for rho != 0
            return np.ones_like(r)
        
        num = r_safe * jb_val**(n_val - 2)
        den = np.sqrt(jb_val**(2 * (n_val - 1)) - 1)
        # We take the absolute value of rho in the calculation for F_rho, 
        # so F_rho is symmetric.
        den_safe = np.where(den == 0, 1e-9, den)
        return num / den_safe

    # Define components of the operator
    rho_mid = (rho[:-1] + rho[1:]) / 2
    
    jb_rho = japanese_bracket(rho)
    jb_rho_mid = japanese_bracket(rho_mid)

    F_rho_val = F_rho_abs(rho, n)
    F_rho_val_mid = F_rho_abs(rho_mid, n)
    
    # Weight function w(rho) for the inner product
    w = jb_rho**(n - 1) * F_rho_val

    # p(rho) coefficient in the differential operator
    p = jb_rho_mid**(n - 1) / F_rho_val_mid

    # Potential V_k(rho)
    V = -k * (k + n - 2) / jb_rho**2 + n * (n - 1) / jb_rho**(2 * n)
    
    # Discretize the operator L_k = (1/w) * d/d_rho(p * d/d_rho) + V
    # We construct the matrix for the generalized eigenvalue problem A*f = lambda*B*f
    # where A is the matrix for p*f''+p'*f' + V*w*f and B is the diagonal matrix for w.
    
    # Matrix for the differential part (1/h^2) * [p_{i-1/2}, -(p_{i-1/2}+p_{i+1/2}), p_{i+1/2}]
    # We impose Dirichlet boundary conditions, f(-rho_max)=f(rho_max)=0, so we operate on internal points.
    internal_points = num_points - 2
    A = np.zeros((internal_points, internal_points))
    
    p_internal = p[1:-1]
    
    A += np.diag(p_internal, k=1)
    A += np.diag(p_internal, k=-1)
    A += np.diag(-(p[1:-1] + p[:-2]))
    A /= h**2
    
    # Add potential term
    V_internal = V[1:-1]
    w_internal = w[1:-1]
    A += np.diag(V_internal * w_internal)
    
    # Create the matrix for the right-hand side of the generalized eigenvalue problem
    B = np.diag(w_internal)

    # Solve the generalized eigenvalue problem
    try:
        eigenvalues = eig(A, B, right=True)[0]
    except np.linalg.LinAlgError:
        print("Eigenvalue computation failed.")
        return 0
    
    # Count positive eigenvalues
    positive_eigenvalues_count = np.sum(eigenvalues.real > 1e-6)
    
    return positive_eigenvalues_count

def solve_problem():
    """
    Solves the problem by checking relevant modes and summing up the number of positive eigenvalues.
    We test for n=3, as the result is expected to be independent of n.
    """
    n_test = 3
    
    # For n=3, V_k is negative for k>=2. So we only need to check k=0 and k=1.
    
    # k=0 mode
    count_k0 = count_positive_eigenvalues(n=n_test, k=0)
    print(f"For n={n_test}, mode k=0 has {count_k0} positive eigenvalue(s).")
    
    # k=1 mode
    count_k1 = count_positive_eigenvalues(n=n_test, k=1)
    print(f"For n={n_test}, mode k=1 has {count_k1} positive eigenvalue(s).")

    # k=2 mode check to confirm analysis
    count_k2 = count_positive_eigenvalues(n=n_test, k=2)
    print(f"For n={n_test}, mode k=2 has {count_k2} positive eigenvalue(s).")

    total_positive_eigenvalues = count_k0 + (n_test-1)*count_k1 # Multiplicity for k=1 modes is (n-1) in S^{n-1}
                                                                # But here we count for the operator L_k once. 
                                                                # All degenerate Y_k have the same radial eqn.
    
    # Total count across all modes is the sum.
    # Spherical harmonics have degeneracies: dim(V_k) > 1 for k>0.
    # We count eigenvalues of L, so we must account for multiplicity.
    # Dim of eigenspace for lambda_k,j is dim(V_k).
    # For S^2 (n=3), dim(V_k) = 2k+1.
    # So total positive eigenvalues are sum over k of (2k+1) * (# pos eig for L_k).
    dim_V0 = 1
    dim_V1 = 2*1+1 # For n=3
    
    total = dim_V0 * count_k0 + dim_V1 * count_k1
    
    print(f"\nThe stability operator L for n={n_test} has a total of {total} positive eigenvalues.")
    # The final equation print
    print("\nThe total number of positive eigenvalues is the sum of counts for each mode, weighted by the degeneracy of the spherical harmonics.")
    print(f"For n={n_test}:")
    print(f"Mode k=0: degeneracy=1, positive eigenvalues={count_k0}")
    print(f"Mode k=1: degeneracy=3, positive eigenvalues={count_k1}")
    print(f"Mode k=2: degeneracy=5, positive eigenvalues={count_k2}")
    print(f"Total = 1 * {count_k0} + 3 * {count_k1} + 5 * {count_k2} + ... = {total}")

solve_problem()