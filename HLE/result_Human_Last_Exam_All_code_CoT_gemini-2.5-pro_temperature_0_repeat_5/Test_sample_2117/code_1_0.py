import numpy as np

def calculate_terms(n):
    """
    Calculates the terms E_P, E_H, S_P, S_H for a given n.
    This function demonstrates that E_P is always zero.
    """
    m = n + 2
    
    # 1. Construct the Cayley-Menger matrix C_n = J - I
    C = np.ones((m, m)) - np.identity(m)

    # 2. Perform the Gaussian Hessenberg Decomposition
    # The algorithm only requires one step for the first column.
    # The transformation T = T_1, where T_1 eliminates entries in column 1.
    # H = T * C * T^{-1}
    # The problem defines C = P * H * P^{-1}, so P = T^{-1}.
    
    T = np.identity(m)
    # Create zeros in the first column (index 0) below the sub-diagonal (index 1)
    # using row 1 as the pivot row.
    for i in range(2, m):
        # T[i, 1] = -C[i, 0] / C[1, 0] = -1 / 1 = -1
        T[i, 1] = -1.0
        
    # P is the inverse of the transformation matrix T
    P = np.linalg.inv(T)
    
    # H is the resulting Hessenberg matrix
    H = T @ C @ np.linalg.inv(T)

    # 3. Calculate E_P
    # Eigenvalues of a unit triangular matrix are all 1.
    eigvals_P = np.linalg.eigvals(P)
    E_P = (np.max(np.real(eigvals_P)) - np.min(np.real(eigvals_P))) / (m - 1)

    # 4. Calculate E_H
    # Eigenvalues of H are the same as C: (n+1) and -1 (n+1 times)
    lambda_max_H = n + 1
    lambda_min_H = -1
    E_H = (lambda_max_H - lambda_min_H) / (m - 1)

    # 5. Calculate S_P
    # S_M = ||M||_F^2 / size
    S_P = np.sum(P**2) / m

    # 6. Calculate S_H
    S_H = np.sum(H**2) / m
    
    return E_P, E_H, S_P, S_H

# Set a sample value for n
n = 3

# Calculate the terms
E_P, E_H, S_P, S_H = calculate_terms(n)
product = E_P * E_H * S_P * S_H

# The symbolic analysis proves the result is always 0.
# We print the equation for n=3 to demonstrate.
print(f"For n = {n}:")
print("The analysis shows that the transformation matrix P is unit triangular, which means all its eigenvalues are 1.")
print(f"Therefore, the average eigenvalue gap of P is E_P = (1 - 1) / ({n+2} - 1) = 0.")
print("\nThe final product is the multiplication of four terms: E_P * E_H * S_P * S_H.")
print("The equation with the calculated values is:")
print(f"{E_P:.4f} * {E_H:.4f} * {S_P:.4f} * {S_H:.4f} = {product:.4f}")
print("\nSince E_P is 0 for all n, the product is always 0.")
print("The least upper bound of the product over all positive integers n is therefore 0.")

final_answer = 0.0
# The final answer is formatted as requested.
# print(f'<<<{final_answer}>>>')