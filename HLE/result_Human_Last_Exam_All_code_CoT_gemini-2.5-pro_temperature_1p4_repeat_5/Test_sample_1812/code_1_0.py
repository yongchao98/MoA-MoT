import numpy as np

def compute_determinant_and_show_steps():
    """
    Computes the determinant of a given 3x3 matrix and prints the calculation steps.
    """
    # The given matrix A
    A = np.array([
        [2, -1, -1],
        [-1, 2, -1],
        [-1, -1, 2]
    ])

    # Unpack the matrix elements for clarity in the formula
    a, b, c = A[0, 0], A[0, 1], A[0, 2]
    d, e, f = A[1, 0], A[1, 1], A[1, 2]
    g, h, i = A[2, 0], A[2, 1], A[2, 2]

    # The determinant formula for a 3x3 matrix via cofactor expansion is:
    # det(A) = a(ei - fh) - b(di - fg) + c(dh - eg)

    # Calculate the sub-determinants (minors of the first row elements)
    sub_det_1 = e * i - f * h
    sub_det_2 = d * i - f * g
    sub_det_3 = d * h - e * g

    # Calculate the full determinant using the formula
    determinant = a * sub_det_1 - b * sub_det_2 + c * sub_det_3

    print("The given matrix A is:")
    print(A)
    
    print("\nThe determinant is calculated using the formula for cofactor expansion along the first row:")
    print("det(A) = a(ei - fh) - b(di - fg) + c(dh - eg)\n")

    print("Step 1: Substitute the values from matrix A into the formula:")
    print(f"det(A) = {a}(({e})*({i}) - ({f})*({h})) - ({b})(({d})*({i}) - ({f})*({g})) + {c}(({d})*({h}) - ({e})*({g}))")

    print("\nStep 2: Calculate the value of each sub-determinant:")
    print(f"det(A) = {a}({sub_det_1}) - ({b})({sub_det_2}) + {c}({sub_det_3})")

    print("\nStep 3: Perform the multiplications:")
    print(f"det(A) = {a * sub_det_1} - ({b * sub_det_2}) + ({c * sub_det_3})")
    
    print("\nStep 4: Compute the final result:")
    print(f"det(A) = {a * sub_det_1 - b * sub_det_2 + c * sub_det_3}")

    print(f"\nThe determinant of the adjacency matrix is: {determinant}")

# Execute the function to print the solution
compute_determinant_and_show_steps()