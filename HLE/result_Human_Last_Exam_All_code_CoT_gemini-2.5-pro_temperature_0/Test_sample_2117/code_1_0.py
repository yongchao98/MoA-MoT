import numpy as np

def solve_for_n(n):
    """
    Calculates the product E_P * E_H * S_P * S_H for a given integer n.
    """
    if not isinstance(n, int) or n <= 0:
        print("Please provide a positive integer for n.")
        return

    print(f"--- Calculation for n = {n} ---")
    m = n + 2

    # 1. Construct the Cayley-Menger matrix C
    C = np.ones((m, m)) - np.identity(m)

    # 2. Perform Gaussian Hessenberg Decomposition C = P * H * P^-1
    # This means H = P^-1 * C * P
    # We find P^-1 by applying elementary transformations to C
    A = C.copy()
    P_inv = np.identity(m)

    for k in range(m - 2):
        # For C = J - I, no permutations are needed.
        # The pivot A[k+1, k] is never zero when an operation is required.
        for i in range(k + 2, m):
            if A[k + 1, k] == 0:
                # This case does not happen for C=J-I when elimination is needed
                continue
            
            multiplier = A[i, k] / A[k + 1, k]
            
            # Create elementary matrix L
            L = np.identity(m)
            L[i, k + 1] = -multiplier
            
            # Apply similarity transformation
            A = L @ A @ np.linalg.inv(L)
            P_inv = L @ P_inv

    # H is the resulting Hessenberg matrix
    H = A
    # P is the inverse of P_inv
    P = np.linalg.inv(P_inv)

    # 3. Calculate E_P, E_H, S_P, S_H
    
    # Eigenvalues of P
    eig_P = np.linalg.eigvals(P)
    # Eigenvalues of H
    eig_H = np.linalg.eigvals(H)

    # E_P: Average eigenvalue gap of P
    # Eigenvalues of P are all 1, so the gap is 0.
    E_P = (np.max(eig_P.real) - np.min(eig_P.real)) / (m - 1)

    # E_H: Average eigenvalue gap of H
    # Eigenvalues of H are n+1 and -1.
    E_H = (np.max(eig_H.real) - np.min(eig_H.real)) / (m - 1)

    # S_P: Mean square of singular values of P
    S_P = np.sum(P**2) / m

    # S_H: Mean square of singular values of H
    S_H = np.sum(H**2) / m

    # 4. Calculate the final product
    product = E_P * E_H * S_P * S_H

    print(f"For n={n}, the matrix size is m={m}.")
    print(f"E_P = (max_eig(P) - min_eig(P)) / (m-1) = ({np.max(eig_P.real):.4f} - {np.min(eig_P.real):.4f}) / {m-1} = {E_P:.4f}")
    print(f"E_H = (max_eig(H) - min_eig(H)) / (m-1) = ({np.max(eig_H.real):.4f} - {np.min(eig_H.real):.4f}) / {m-1} = {E_H:.4f}")
    print(f"S_P = ||P||_F^2 / m = {np.sum(P**2):.4f} / {m} = {S_P:.4f}")
    print(f"S_H = ||H||_F^2 / m = {np.sum(H**2):.4f} / {m} = {S_H:.4f}")
    print("\nFinal Equation:")
    print(f"E_P * E_H * S_P * S_H = {E_P:.4f} * {E_H:.4f} * {S_P:.4f} * {S_H:.4f} = {product:.4f}")
    print("-" * 20)


# Run the calculation for a sample value, e.g., n=3
solve_for_n(3)

print("\nAs demonstrated, E_P is always 0 because the Gaussian Hessenberg decomposition")
print("of the given Cayley-Menger matrix results in a unit lower triangular matrix P,")
print("whose eigenvalues are all 1. This makes the entire product zero for any n.")
print("The set of values for the product over all positive integers n is {0}.")
print("Therefore, the least upper bound of this set is 0.")
