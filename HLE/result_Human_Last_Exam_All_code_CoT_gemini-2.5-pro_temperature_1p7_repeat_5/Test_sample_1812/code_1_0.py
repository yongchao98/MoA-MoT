import numpy as np

def compute_determinant():
    """
    Computes the determinant of the given 3x3 adjacency matrix
    and prints the step-by-step calculation.
    """
    # The matrix A is provided in the problem description.
    # The mention of the Markov quiver and m_{1/2}=13 is context;
    # the core task is to find the determinant of the given matrix A.
    A = np.array([
        [2, -1, -1],
        [-1, 2, -1],
        [-1, -1, 2]
    ])

    print("The given adjacency matrix is:")
    print(A)
    print("-" * 30)

    # Extract elements for printing the equation clearly
    a, b, c = A[0, 0], A[0, 1], A[0, 2]
    d, e, f = A[1, 0], A[1, 1], A[1, 2]
    g, h, i = A[2, 0], A[2, 1], A[2, 2]

    # The formula for the determinant of a 3x3 matrix is:
    # det(A) = a(ei - fh) - b(di - fg) + c(dh - eg)

    print("To compute the determinant, we substitute the matrix elements into the formula.")
    print("\nFinal Equation:")
    # Print the full equation with all numbers from the matrix
    print(f"det(A) = ({a}) * (({e})*({i}) - ({f})*({h})) - ({b}) * (({d})*({i}) - ({f})*({g})) + ({c}) * (({d})*({h}) - ({e})*({g}))")

    # Step-by-step calculation
    term1_inner = e * i - f * h
    term2_inner = d * i - f * g
    term3_inner = d * h - e * g
    print("\nStep 1: Calculate the terms in parentheses")
    print(f"det(A) = ({a}) * ({term1_inner}) - ({b}) * ({term2_inner}) + ({c}) * ({term3_inner})")

    term1 = a * term1_inner
    term2 = -b * term2_inner
    term3 = c * term3_inner
    print("\nStep 2: Perform the multiplications")
    print(f"det(A) = ({term1}) + ({term2}) + ({term3})")

    # Final result from the step-by-step calculation
    final_determinant = term1 + term2 + term3
    print("\nStep 3: Sum the terms to get the final determinant")
    print(f"det(A) = {final_determinant}")
    print("-" * 30)

    # Verify with numpy's determinant function
    det_A_numpy = np.linalg.det(A)
    print(f"The determinant as calculated by numpy is: {det_A_numpy:.0f}")

if __name__ == "__main__":
    compute_determinant()