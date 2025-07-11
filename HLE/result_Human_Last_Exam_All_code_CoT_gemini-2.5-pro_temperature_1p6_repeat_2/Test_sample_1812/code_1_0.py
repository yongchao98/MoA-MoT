import numpy as np

def solve():
    """
    Computes the determinant of the adjacency matrix of the Markov quiver.
    """
    # The matrix is given as:
    # A = [[2, -1, -1],
    #      [-1, 2, -1],
    #      [-1, -1, 2]]
    # The reference to m_{1/2} = 13 is context for the origin of the matrix,
    # and does not change the matrix itself.
    A = np.array([
        [2, -1, -1],
        [-1, 2, -1],
        [-1, -1, 2]
    ])

    print("The given matrix A is:")
    print(A)
    print("\nThe determinant of a 3x3 matrix [[a, b, c], [d, e, f], [g, h, i]] is calculated using the formula:")
    print("det(A) = a(ei - fh) - b(di - fg) + c(dh - eg)")

    a, b, c = A[0, 0], A[0, 1], A[0, 2]
    d, e, f = A[1, 0], A[1, 1], A[1, 2]
    g, h, i = A[2, 0], A[2, 1], A[2, 2]

    # Calculate determinant using the formula to show the steps
    term1_val = e * i - f * h
    term2_val = d * i - f * g
    term3_val = d * h - e * g

    det_val = a * term1_val - b * term2_val + c * term3_val

    print("\nSubstituting the values from matrix A into the formula:")
    print(f"det(A) = ({a}) * (({e})*({i}) - ({f})*({h})) - ({b}) * (({d})*({i}) - ({f})*({g})) + ({c}) * (({d})*({h}) - ({e})*({g}))")
    print(f"det(A) = ({a}) * ({term1_val}) - ({b}) * ({term2_val}) + ({c}) * ({term3_val})")
    print(f"det(A) = ({a * term1_val}) - ({b * term2_val}) + ({c * term3_val})")
    print(f"det(A) = {a * term1_val - b * term2_val + c * term3_val}")
    
    # Also compute with numpy for verification
    det_A_numpy = np.linalg.det(A)
    
    print(f"\nThe final computed determinant is: {int(round(det_A_numpy))}")

solve()