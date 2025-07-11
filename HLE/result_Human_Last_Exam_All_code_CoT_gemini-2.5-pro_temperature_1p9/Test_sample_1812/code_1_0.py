import numpy as np

def compute_determinant():
    """
    Computes and prints the step-by-step determinant calculation for the given 3x3 matrix.
    """
    # The adjacency matrix of the Markov quiver is given as A.
    A = np.array([
        [2, -1, -1],
        [-1, 2, -1],
        [-1, -1, 2]
    ])

    # For the formula, we can represent the matrix as:
    # A = [[a, b, c],
    #      [d, e, f],
    #      [g, h, i]]
    a, b, c = A[0, 0], A[0, 1], A[0, 2]
    d, e, f = A[1, 0], A[1, 1], A[1, 2]
    g, h, i = A[2, 0], A[2, 1], A[2, 2]

    # The determinant formula for a 3x3 matrix is:
    # det(A) = a(ei - fh) - b(di - fg) + c(dh - eg)
    
    print("The matrix A is:")
    print(A)
    print("\nTo compute the determinant, we use the formula:")
    print("det(A) = a(ei - fh) - b(di - fg) + c(dh - eg)")

    print("\nSubstituting the values from matrix A into the formula:")
    final_equation_str = (
        f"det(A) = ({a}) * (({e})*({i}) - ({f})*({h}))"
        f" - ({b}) * (({d})*({i}) - ({f})*({g}))"
        f" + ({c}) * (({d})*({h}) - ({e})*({g}))"
    )
    print(final_equation_str)

    # Perform the calculations step by step
    val_ei_fh = e * i - f * h
    val_di_fg = d * i - f * g
    val_dh_eg = d * h - e * g
    
    print(f"\n       = ({a})*({val_ei_fh}) - ({b})*({val_di_fg}) + ({c})*({val_dh_eg})")
    
    term1 = a * val_ei_fh
    term2 = -b * val_di_fg
    term3 = c * val_dh_eg
    
    print(f"       = {term1} + {term2} + {term3}")

    determinant = np.linalg.det(A)
    print(f"       = {int(determinant)}")
    print("\nTherefore, the determinant of the given matrix is 0.")

# Execute the function to show the result
compute_determinant()