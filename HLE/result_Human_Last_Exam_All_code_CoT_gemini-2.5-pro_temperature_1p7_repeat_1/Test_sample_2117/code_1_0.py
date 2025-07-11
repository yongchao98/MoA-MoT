import numpy as np

def solve(n):
    """
    Calculates the product E_P * E_H * S_P * S_H for a given n.
    """
    if n <= 0:
        print("n must be a positive integer.")
        return

    # 1. Construct the Cayley-Menger matrix C
    size = n + 1
    C = np.ones((size, size)) - np.identity(size)
    
    # Keep the original C for verification
    C_orig = C.copy()

    # 2. Perform Gaussian Hessenberg Decomposition: H = M C M^-1
    # We will find P such that C = P H P^-1, so P = M^-1
    # Start with P = I, and accumulate the inverse transformations.
    M_inv_total = np.identity(size)
    H = C.copy()

    for k in range(size - 2):
        for i in range(k + 2, size):
            # The pivot is H[k+1, k]. For C=J-I, it can be shown
            # that pivots will be non-zero.
            if H[k + 1, k] == 0:
                # This case should not be reached for J-I matrix
                raise ValueError("Pivot is zero, pivoting is required.")

            multiplier = H[i, k] / H[k + 1, k]
            
            # Create elementary matrix M and its inverse M_inv
            M = np.identity(size)
            M[i, k + 1] = -multiplier
            M_inv = np.identity(size)
            M_inv[i, k + 1] = multiplier
            
            # Update H with the similarity transformation
            H = M @ H @ M_inv
            
            # Accumulate the total transformation P = M_0^-1 * M_1^-1 * ...
            M_inv_total = M_inv_total @ M_inv

    P = M_inv_total
    
    # 3. Calculate all four quantities

    # For H
    eig_H, _ = np.linalg.eig(H)
    eig_H_real = np.real(eig_H)
    E_H = (np.max(eig_H_real) - np.min(eig_H_real)) / (size - 1)
    S_H = np.sum(np.linalg.svd(H, compute_uv=False)**2) / size

    # For P
    eig_P, _ = np.linalg.eig(P)
    eig_P_real = np.real(eig_P)
    # The crucial step: P is unitriangular, so its eigenvalues are all 1.
    E_P = (np.max(eig_P_real) - np.min(eig_P_real)) / (size - 1)
    S_P = np.sum(np.linalg.svd(P, compute_uv=False)**2) / size

    product = E_P * E_H * S_P * S_H

    print(f"For n = {n}:")
    print(f"E_P = ({np.max(eig_P_real):.4f} - {np.min(eig_P_real):.4f}) / {size - 1} = {E_P:.4f}")
    print(f"E_H = ({np.max(eig_H_real):.4f} - {np.min(eig_H_real):.4f}) / {size - 1} = {E_H:.4f}")
    print(f"S_P = {S_P:.4f}")
    print(f"S_H = {S_H:.4f}")
    print("\nFinal Product Calculation:")
    print(f"E_P * E_H * S_P * S_H = {E_P:.4f} * {E_H:.4f} * {S_P:.4f} * {S_H:.4f} = {product:.4f}")
    
# Let's run for n=3 as an example
solve(3)

print("\nConclusion:")
print("The analysis shows E_P is always 0 for any positive integer n.")
print("Therefore, the product E_P * E_H * S_P * S_H is always 0.")
print("The least upper bound of a set of zeros is 0.")
<<<0>>>