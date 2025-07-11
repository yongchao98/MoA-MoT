import numpy as np

def calculate_metrics(matrix):
    """Calculates the average eigenvalue gap and mean square of singular values."""
    # Ensure matrix is real for eigenvalue sorting
    if np.iscomplexobj(matrix):
        eigenvalues = np.linalg.eigvals(matrix)
        real_parts = np.real(eigenvalues)
        real_parts.sort()
        # Using real parts for eigenvalue gap definition with complex numbers
        avg_eigenvalue_gap = (real_parts[-1] - real_parts[0]) / (len(real_parts) - 1)
    else:
        eigenvalues = np.linalg.eigvals(matrix)
        eigenvalues.sort()
        avg_eigenvalue_gap = (eigenvalues[-1] - eigenvalues[0]) / (len(eigenvalues) - 1)
    
    singular_values = np.linalg.svd(matrix, compute_uv=False)
    mean_square_singular_values = np.mean(singular_values**2)
    
    return avg_eigenvalue_gap, mean_square_singular_values

def main():
    """
    Solves the problem for n=1 to find the least upper bound.
    """
    n = 1
    k = n + 2

    # 1. Define the Cayley-Menger matrix M for n=1
    M = np.ones((k, k)) - np.identity(k)

    # 2. Construct a symmetric orthogonal matrix P of eigenvectors of M
    # The first eigenvector for eigenvalue (n+1)=2 is (1,1,1)/sqrt(3)
    # We construct an orthonormal basis for the eigenspace for eigenvalue -1
    # such that the resulting matrix P is symmetric.
    # P_ij = u_j(i) where u_j is the j-th eigenvector (column).
    # For P to be symmetric, u_j(i) = u_i(j).
    # This leads to a specific choice of eigenvectors.
    d = (3 - np.sqrt(3)) / 6
    e = (-3 - np.sqrt(3)) / 6
    
    P = np.array([
        [1/np.sqrt(3), 1/np.sqrt(3), 1/np.sqrt(3)],
        [1/np.sqrt(3), d, e],
        [1/np.sqrt(3), e, d]
    ])

    # 3. Define the Hessenberg matrix H
    # Since P is orthogonal, P_inv = P.T
    # H is diagonal because P consists of eigenvectors of M.
    H = P.T @ M @ P

    # 4. Calculate the four metrics
    E_P, S_P = calculate_metrics(P)
    E_H, S_H = calculate_metrics(H)
    
    # 5. Calculate the final product
    product = E_P * E_H * S_P * S_H

    # The trace of this P is 1, and since it's symmetric and orthogonal,
    # its eigenvalues must be {1, 1, -1}.
    # E_P = (1 - (-1)) / (3-1) = 1
    # S_P = 1 since it's orthogonal
    # Eigenvalues of H are {2, -1, -1}
    # E_H = (2 - (-1)) / (3-1) = 1.5
    # S_H = (2^2 + (-1)^2 + (-1)^2) / 3 = (4+1+1)/3 = 2
    # Product = 1 * 1.5 * 1 * 2 = 3
    
    print("For n=1, we can construct a specific decomposition that yields the maximum value.")
    print("This demonstrates that the least upper bound is achievable.")
    print("\nCalculating the components for the product E_P * E_H * S_P * S_H:")
    print(f"E_P (Average Eigenvalue Gap of P) = {E_P:.4f}")
    print(f"E_H (Average Eigenvalue Gap of H) = {E_H:.4f}")
    print(f"S_P (Mean Square of Singular Values of P) = {S_P:.4f}")
    print(f"S_H (Mean Square of Singular Values of H) = {S_H:.4f}")
    print("\nThe final equation is:")
    print(f"{E_P:.4f} * {E_H:.4f} * {S_P:.4f} * {S_H:.4f} = {product:.4f}")
    
    print("\nBased on the theoretical analysis and this numerical confirmation, the least upper bound is 3.")

if __name__ == "__main__":
    main()
