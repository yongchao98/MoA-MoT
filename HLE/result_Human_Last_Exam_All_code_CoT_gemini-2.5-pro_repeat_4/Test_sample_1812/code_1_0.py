import numpy as np

def solve_determinant():
    """
    Computes the determinant of the adjacency matrix of the Markov quiver.
    """
    # The given matrix A
    A = np.array([
        [2, -1, -1],
        [-1, 2, -1],
        [-1, -1, 2]
    ])

    # Extract elements for the determinant formula
    a, b, c = A[0, 0], A[0, 1], A[0, 2]
    d, e, f = A[1, 0], A[1, 1], A[1, 2]
    g, h, i = A[2, 0], A[2, 1], A[2, 2]

    print("The given matrix is:")
    print(A)
    print("\nTo compute the determinant, we use the expansion formula for a 3x3 matrix:")
    print("det(A) = a(ei - fh) - b(di - fg) + c(dh - eg)")

    print("\nSubstituting the values from the matrix into the formula:")
    # Print the equation with all numbers
    print(f"det(A) = ({a}) * (({e})*({i}) - ({f})*({h})) - ({b}) * (({d})*({i}) - ({f})*({g})) + ({c}) * (({d})*({h}) - ({e})*({g}))")

    # Calculate intermediate values within the parentheses
    sub_term1 = e * i - f * h
    sub_term2 = d * i - f * g
    sub_term3 = d * h - e * g
    print(f"       = {a} * ({sub_term1}) - ({b}) * ({sub_term2}) + ({c}) * ({sub_term3})")

    # Calculate the products
    term1 = a * sub_term1
    term2 = -b * sub_term2
    term3 = c * sub_term3
    print(f"       = {term1} + ({term2}) + ({term3})")

    # Calculate the final determinant
    determinant = term1 + term2 + term3
    print(f"       = {determinant}")

    # Verify with numpy's built-in function
    final_determinant = np.linalg.det(A)
    print(f"\nThe final computed determinant is {int(round(final_determinant))}.")

solve_determinant()
<<<0>>>