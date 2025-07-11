import numpy as np

def solve():
    """
    This function solves the problem by demonstrating that the product is always 0.
    We will choose a sample value for n, e.g., n=4, and compute all the quantities.
    The logic holds for any positive integer n.
    """
    # The user can change n to any positive integer to test.
    n = 4
    m = n + 2  # Size of the Cayley-Menger matrix

    # Step 1: Construct the Cayley-Menger matrix C.
    C = np.ones((m, m)) - np.identity(m)

    # Step 2: Define the transformation matrix P from the Gaussian Hessenberg decomposition.
    # For the specific structure of C, the decomposition algorithm produces a
    # unit lower triangular matrix P with 1s on the diagonal and 1s in the
    # second column below the diagonal.
    P = np.identity(m)
    P[2:, 1] = 1

    # The inverse of P is needed to find H.
    P_inv = np.identity(m)
    P_inv[2:, 1] = -1

    # The Hessenberg matrix H is calculated as H = P_inv * C * P.
    H = P_inv @ C @ P

    # Step 3: Calculate the average eigenvalue gap E_M for P and H.
    # E_M = (lambda_max - lambda_min) / (k-1)
    
    # For P:
    eigvals_P = np.linalg.eigvals(P)
    # The eigenvalues of a triangular matrix are its diagonal entries, which are all 1.
    # So, max and min eigenvalues are both 1.
    E_P = (np.max(eigvals_P.real) - np.min(eigvals_P.real)) / (m - 1)

    # For H: H is similar to C, so they have the same eigenvalues.
    eigvals_H = np.linalg.eigvals(H)
    E_H = (np.max(eigvals_H.real) - np.min(eigvals_H.real)) / (m - 1)
    
    # Step 4: Calculate the mean square of singular values S_M for P and H.
    # S_M = (1/k) * ||M||_F^2 = (1/k) * sum(abs(elements)^2)
    
    # For P:
    S_P = np.sum(np.square(np.abs(P))) / m

    # For H:
    S_H = np.sum(np.square(np.abs(H))) / m
    
    # Step 5: Calculate the final product
    product = E_P * E_H * S_P * S_H

    print(f"This script calculates the product E_P * E_H * S_P * S_H for n={n}.")
    print("----------------------------------------------------------")
    print(f"The final equation is E_P * E_H * S_P * S_H")
    print(f"E_P = {E_P:.6f}")
    print(f"E_H = {E_H:.6f}")
    print(f"S_P = {S_P:.6f}")
    print(f"S_H = {S_H:.6f}")
    print(f"Product = {E_P:.6f} * {E_H:.6f} * {S_P:.6f} * {S_H:.6f} = {product:.6f}")
    print("----------------------------------------------------------")
    print("\nAs shown by the calculation, E_P is 0. This holds true for any n.")
    print("Because E_P is always 0, the entire product is always 0.")
    print("The least upper bound of a set containing only 0 is 0.")

solve()
<<<0>>>