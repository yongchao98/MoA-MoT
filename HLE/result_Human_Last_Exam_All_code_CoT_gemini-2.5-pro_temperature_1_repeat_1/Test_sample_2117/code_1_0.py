import numpy as np

def solve_for_n(n):
    """
    Calculates the terms E_P, E_H, S_P, S_H for a given integer n
    and prints the final product equation.
    """
    # The size of the matrices is m x m where m = n + 2
    m = n + 2

    # Step 1: Determine the matrices P and H based on the decomposition
    
    # P is a unit lower triangular matrix derived from the Gauss transforms.
    # P = I + sum_{i=3 to m} E_{i,2}
    P = np.identity(m)
    if m > 2:
        P[2:, 1] = 1

    # H is the upper Hessenberg matrix from the decomposition.
    # H = P^{-1} * M_n * P, where M_n is the Cayley-Menger matrix.
    H = np.zeros((m, m))
    if m > 0:
        H[0, 1] = n + 1
        H[0, 2:] = 1
    if m > 1:
        H[1, 0] = 1
        H[1, 1] = n
        H[1, 2:] = 1
    if m > 2:
        for i in range(2, m):
            H[i, i] = -1

    # Step 2: Calculate E_P and S_P for matrix P

    # Eigenvalues of a triangular matrix are its diagonal entries.
    # All diagonal entries of P are 1.
    eigvals_P = np.linalg.eigvals(P)
    eigvals_P_sorted = np.sort(np.real(eigvals_P))
    
    # Average eigenvalue gap for P
    e_p = np.mean(np.abs(np.diff(eigvals_P_sorted)))
    
    # Mean square of singular values for P
    s_p = np.mean(np.linalg.svd(P, compute_uv=False)**2)

    # Step 3: Calculate E_H and S_H for matrix H

    # Eigenvalues of H are n+1 (multiplicity 1) and -1 (multiplicity n+1)
    eigvals_H = np.full(m, -1.0)
    eigvals_H[0] = n + 1.0
    eigvals_H_sorted = np.sort(eigvals_H)

    # Average eigenvalue gap for H
    # The only non-zero gap is between -1 and n+1, which is (n+1) - (-1) = n+2.
    # The sum of gaps is n+2, and there are m-1 = n+1 gaps.
    e_h = (n + 2.0) / (n + 1.0)

    # Mean square of singular values for H
    # S_H = ||H||_F^2 / m = (2n^2 + 5n + 2) / (n+2) = 2n + 1
    s_h = np.mean(np.linalg.svd(H, compute_uv=False)**2)
    
    # Step 4: Calculate the final product
    
    product = e_p * e_h * s_p * s_h

    # Print the results in the required format
    print(f"For n = {n}:")
    print("The final equation is E_P * E_H * S_P * S_H = Product")
    print(f"Each number in the final equation:")
    print(f"{e_p} * {e_h} * {s_p} * {s_h} = {product}")
    print("-" * 20)

def main():
    """
    Main function to explain the result and run the calculation.
    """
    print("Step-by-step derivation showed that for any n, the matrix P is unit lower triangular.")
    print("This means all eigenvalues of P are 1.")
    print("Consequently, the average eigenvalue gap of P, denoted E_P, is always 0.")
    print("Therefore, the entire product E_P * E_H * S_P * S_H is always 0.")
    print("The least upper bound of a sequence of zeros is 0.\n")
    print("Let's demonstrate this with a numerical example for n=4.")
    
    # Run the calculation for a sample n
    solve_for_n(4)

    print("As shown, E_P is 0, which makes the product 0.")
    print("This holds true for any positive integer n.")
    print("Thus, the least upper bound of the product is 0.")

if __name__ == "__main__":
    main()