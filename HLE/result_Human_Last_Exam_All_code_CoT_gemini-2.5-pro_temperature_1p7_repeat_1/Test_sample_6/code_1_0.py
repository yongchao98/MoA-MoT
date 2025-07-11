import numpy as np
from scipy.linalg import eigh

def count_kk_modes():
    """
    Calculates the number of Kaluza-Klein graviton modes with mass squared below 14.
    """
    # 1. Discretize the domain [0, 2*pi]
    N = 2000  # Number of grid points for accuracy
    x = np.linspace(0, 2 * np.pi, N, endpoint=False)
    h = 2 * np.pi / N  # Grid spacing

    # 2. Define the warp factor and the coefficient functions p(x) and w(x)
    def A(x_val):
        """ The warp factor function A(x) """
        return np.sin(x_val) + 4 * np.cos(x_val)
    
    def p(x_val):
        """ The coefficient function p(x) in the Sturm-Liouville operator """
        return np.exp(3 * A(x_val))
    
    def w(x_val):
        """ The weight function w(x) in the Sturm-Liouville problem """
        return np.exp(5 * A(x_val))

    # 3. Construct the finite difference matrices K and M
    # The equation is K*psi = lambda*M*psi, where lambda = m^2.
    
    # Values of p(x) at the midpoints of the grid intervals
    p_half = p(x + h / 2)
    
    # Assemble the stiffness matrix K using centered differences
    # Diagonal terms
    K_diag = (p_half + np.roll(p_half, 1)) / h**2
    # Off-diagonal terms
    K_offdiag = -p_half / h**2
    
    # Construct K as a circulant matrix for periodic boundary conditions
    K = np.diag(K_diag) + np.diag(K_offdiag[:-1], k=1) + np.diag(np.roll(K_offdiag, 1)[1:], k=-1)
    K[0, N-1] = K_offdiag[-1] # Correct corner element for K_0,N-1
    K[N-1, 0] = K_offdiag[-1] # Correct corner element for K_N-1,0

    # Assemble the mass matrix M (which is diagonal)
    M = np.diag(w(x))

    # 4. Solve the generalized eigenvalue problem
    # eigh is used for Hermitian matrices; K and M are real and symmetric.
    eigenvalues = eigh(K, M, eigvals_only=True)
    eigenvalues = np.sort(eigenvalues)

    # 5. Count the number of eigenvalues below the threshold of 14
    threshold = 14
    eigenvalues_below_threshold = eigenvalues[eigenvalues < threshold]
    count = len(eigenvalues_below_threshold)

    # Output the results
    print("The mass eigenvalues (m^2) are determined by the Sturm-Liouville equation:")
    print("-d/dx( p(x) * dPsi/dx ) = m^2 * w(x) * Psi")
    print("with p(x) = exp(3*A(x)), w(x) = exp(5*A(x)), and A(x) = sin(x) + 4*cos(x).")
    print("\nWe solve this numerically. The calculated eigenvalues m^2 below 14 are:")
    # Grouping nearly identical values to show degeneracy
    groups = []
    if count > 0:
        current_group = [eigenvalues_below_threshold[0]]
        for i in range(1, count):
            if np.isclose(eigenvalues_below_threshold[i], current_group[0], rtol=1e-5):
                current_group.append(eigenvalues_below_threshold[i])
            else:
                groups.append(current_group)
                current_group = [eigenvalues_below_threshold[i]]
        groups.append(current_group)
    
    print_groups = []
    sum_str_parts = []
    for group in groups:
        print_groups.append(f"m^2 ~ {group[0]:.4f} (degeneracy {len(group)})")
        sum_str_parts.append(str(len(group)))
        
    print("\n".join(print_groups))
    
    final_equation = " + ".join(sum_str_parts)
    print(f"\nThe total count is the sum of these degeneracies:")
    print(f"Total Count = {final_equation} = {count}")


count_kk_modes()

# The final answer as an integer value
<<<13>>>