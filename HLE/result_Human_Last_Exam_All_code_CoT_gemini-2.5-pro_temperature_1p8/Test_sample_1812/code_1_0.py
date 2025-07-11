import numpy as np

def compute_determinant():
    """
    Computes the determinant of the given adjacency matrix and prints the calculation.
    """
    # Define the adjacency matrix A
    A = np.array([
        [2, -1, -1],
        [-1, 2, -1],
        [-1, -1, 2]
    ])

    print("The adjacency matrix is:")
    print(A)

    # Unpack the matrix elements for displaying the calculation
    a, b, c = A[0, 0], A[0, 1], A[0, 2]
    d, e, f = A[1, 0], A[1, 1], A[1, 2]
    g, h, i = A[2, 0], A[2, 1], A[2, 2]
    
    # Calculate the determinant using numpy
    determinant = np.linalg.det(A)

    # Print the final equation with each number from the matrix
    # The determinant formula is: a(ei - fh) - b(di - fg) + c(dh - eg)
    print("\nThe determinant calculation is shown by the following equation:")
    print(f"{a} * (({e}) * ({i}) - ({f}) * ({h})) - ({b}) * (({d}) * ({i}) - ({f}) * ({g})) + ({c}) * (({d}) * ({h}) - ({e}) * ({g}))")
    
    # Calculate and show the intermediate steps
    term1 = e * i - f * h
    term2 = d * i - f * g
    term3 = d * h - e * g
    print(f"= {a} * ({term1}) - ({b}) * ({term2}) + ({c}) * ({term3})")
    
    # Show the final step of the arithmetic
    print(f"= {a * term1} + {b * term2 * -1} + {c * term3}")

    print("\nFinal determinant:")
    # np.linalg.det returns a float, so we round it and cast to int for a clean answer
    print(int(round(determinant)))

if __name__ == '__main__':
    compute_determinant()
<<<0>>>