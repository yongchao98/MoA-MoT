import numpy as np

def solve():
    """
    Computes the determinant of the adjacency matrix of the Markov quiver.
    """
    # The matrix A is given as:
    A = np.array([
        [2, -1, -1],
        [-1, 2, -1],
        [-1, -1, 2]
    ])

    print("The adjacency matrix A is:")
    print(A)
    print("")

    # For a 3x3 matrix:
    # [[a, b, c],
    #  [d, e, f],
    #  [g, h, i]]
    # The determinant is a(ei - fh) - b(di - fg) + c(dh - eg)
    a, b, c = A[0, 0], A[0, 1], A[0, 2]
    d, e, f = A[1, 0], A[1, 1], A[1, 2]
    g, h, i = A[2, 0], A[2, 1], A[2, 2]

    # Calculate the determinant
    determinant = np.linalg.det(A)

    # To satisfy the output requirement, we show the calculation.
    print("The determinant calculation is as follows:")
    
    # Step 1: Show the formula with substituted values
    calc_str1 = f"det(A) = ({a}) * (({e})*({i}) - ({f})*({h})) - ({b}) * (({d})*({i}) - ({f})*({g})) + ({c}) * (({d})*({h}) - ({e})*({g}))"
    print(calc_str1)

    # Step 2: Show the result of the cofactor calculations
    term1 = e*i - f*h
    term2 = d*i - f*g
    term3 = d*h - e*g
    calc_str2 = f"det(A) = ({a}) * ({term1}) - ({b}) * ({term2}) + ({c}) * ({term3})"
    print(calc_str2)

    # Step 3: Show the multiplied terms
    final_term1 = a * term1
    final_term2 = -b * term2
    final_term3 = c * term3
    calc_str3 = f"det(A) = {final_term1} + ({final_term2}) + ({final_term3})"
    print(calc_str3)

    # Step 4: Show the final result
    print(f"det(A) = {final_term1 + final_term2 + final_term3}")
    
    # Using integer for clean output as determinant is an integer.
    print(f"\nThe determinant of the matrix is: {int(determinant)}")

solve()