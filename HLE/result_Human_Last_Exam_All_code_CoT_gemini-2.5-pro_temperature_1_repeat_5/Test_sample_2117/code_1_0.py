import numpy as np

def E_M(M):
    """
    Computes the average eigenvalue gap of a matrix M.
    For matrices with complex eigenvalues (like orthogonal matrices), 
    we sort by the real part as a standard convention.
    """
    eigenvalues = np.sort(np.real(np.linalg.eigvals(M)))
    m = len(eigenvalues)
    if m < 2:
        return 0
    return (eigenvalues[-1] - eigenvalues[0]) / (m - 1)

def S_M(M):
    """
    Computes the mean square of the singular values of M.
    For a real matrix M, this is (1/m) * ||M||_F^2.
    """
    m = M.shape[0]
    # For an orthogonal matrix, ||M||_F^2 = m, so S_M = 1.
    # For a symmetric matrix, ||M||_F^2 = Tr(M^T M) = Tr(M^2) = sum(lambda_i^2).
    # S_M = (1/m) * sum(lambda_i^2)
    # Using np.square on the matrix is equivalent to ||M||_F^2 for real matrices.
    return (1 / m) * np.sum(np.square(M))

def solve():
    """
    Calculates the least upper bound by evaluating the case for n=1,
    which maximizes the derived expression 2 + 2/n.
    """
    # The least upper bound is achieved for n=1.
    n = 1
    m = n + 1  # Matrix dimension

    # Define the Cayley-Menger matrix for n=1
    C_M = np.ones((m, m)) - np.identity(m)

    # We use the spectral decomposition C_M = P * H * P^-1.
    # H is the diagonal matrix of eigenvalues (which is a valid Hessenberg matrix).
    # P is the orthogonal matrix of eigenvectors.
    eigenvalues, P = np.linalg.eigh(C_M)
    H = np.diag(eigenvalues)

    # Calculate the four quantities for this n=1 decomposition
    E_P = E_M(P)
    S_P = S_M(P)
    E_H = E_M(H)
    S_H = S_M(H)

    # Calculate the final product
    product = E_P * E_H * S_P * S_H

    print("The least upper bound for the product is achieved when n=1.")
    print(f"For n=1, the equation E_P * E_H * S_P * S_H = Product gives:")
    print(f"{E_P:.4f} * {E_H:.4f} * {S_P:.4f} * {S_H:.4f} = {product:.4f}")

solve()
<<<4>>>