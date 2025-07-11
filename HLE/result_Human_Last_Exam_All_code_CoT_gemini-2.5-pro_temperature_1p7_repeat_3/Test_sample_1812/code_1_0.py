import numpy as np

def compute_determinant():
    """
    Computes the determinant of the given 3x3 matrix and shows the calculation steps.
    """
    # The adjacency matrix of the Markov quiver is given as:
    # A = [[2, -1, -1], [-1, 2, -1], [-1, -1, 2]]
    # Note: The context of the "Markov quiver for m_1/2 = 13" is extra information.
    # The core task is to compute the determinant of the explicitly provided matrix A.
    
    A = np.array([
        [2, -1, -1],
        [-1, 2, -1],
        [-1, -1, 2]
    ])

    # For a 3x3 matrix [[a, b, c], [d, e, f], [g, h, i]], the determinant is
    # calculated using the formula: a(ei - fh) - b(di - fg) + c(dh - eg).
    
    # Extract elements from the matrix
    a, b, c = A[0, 0], A[0, 1], A[0, 2]
    d, e, f = A[1, 0], A[1, 1], A[1, 2]
    g, h, i = A[2, 0], A[2, 1], A[2, 2]

    # Calculate each term of the expansion
    term1 = a * (e * i - f * h)
    term2 = -b * (d * i - f * g)
    term3 = c * (d * h - e * g)
    
    # Sum the terms to get the final determinant
    determinant = term1 + term2 + term3

    print("The matrix A is:")
    print(A)
    print("\nThe determinant is computed using the cofactor expansion formula:")
    print("det(A) = a(ei - fh) - b(di - fg) + c(dh - eg)")
    
    print("\nStep 1: Substitute the values from matrix A into the formula.")
    print(f"det(A) = {a}(({e})*({i}) - ({f})*({h})) - ({b})(({d})*({i}) - ({f})*({g})) + ({c})(({d})*({h}) - ({e})*({g}))")
    
    print("\nStep 2: Calculate the values inside the parentheses (the 2x2 sub-determinants).")
    print(f"det(A) = {a}({e*i - f*h}) - ({b})({d*i - f*g}) + ({c})({d*h - e*g})")

    print("\nStep 3: Perform the multiplications.")
    print(f"det(A) = ({term1}) + ({term2}) + ({term3})")

    print("\nStep 4: Calculate the final sum.")
    print(f"det(A) = {determinant}")

compute_determinant()