import numpy as np

def solve_determinant():
    """
    Computes and explains the determinant of the given 3x3 matrix.
    """
    # The problem mentions a Markov quiver for m_{1/2} = 13, but the core task is to
    # compute the determinant of the explicitly provided matrix A.
    A = np.array([
        [2, -1, -1],
        [-1, 2, -1],
        [-1, -1, 2]
    ])

    print("The given adjacency matrix is:")
    print(A)
    print("\nThe formula for the determinant of a 3x3 matrix [[a, b, c], [d, e, f], [g, h, i]] is:")
    print("det = a(ei - fh) - b(di - fg) + c(dh - eg)")
    
    # Extract elements for printing the formula
    a, b, c = A[0, 0], A[0, 1], A[0, 2]
    d, e, f = A[1, 0], A[1, 1], A[1, 2]
    g, h, i = A[2, 0], A[2, 1], A[2, 2]

    print("\nSubstituting the values from matrix A into the formula:")
    # Print the equation with numbers
    print(f"det(A) = {a}(({e} * {i}) - ({f} * {h})) - ({b})(({d} * {i}) - ({f} * {g})) + ({c})(({d} * {h}) - ({e} * {g}))")

    # Calculate intermediate terms
    term1 = e * i - f * h
    term2 = d * i - f * g
    term3 = d * h - e * g
    
    print("\nSimplifying the terms inside the parentheses:")
    print(f"det(A) = {a}({term1}) - ({b})({term2}) + ({c})({term3})")

    # Calculate the value of each major term
    term1_val = a * term1
    term2_val = b * term2
    term3_val = c * term3
    
    print("\nMultiplying the terms:")
    print(f"det(A) = {term1_val} - ({term2_val}) + ({term3_val})")

    # Final calculation
    final_result = term1_val - term2_val + term3_val
    print("\nFinal calculation:")
    print(f"det(A) = {final_result}")

    # Using numpy to verify the result
    det_A = np.linalg.det(A)
    # Round to handle potential floating point inaccuracies for what should be an integer result
    final_det = int(np.round(det_A))
    
    print(f"\nTherefore, the determinant of the adjacency matrix is {final_det}.")

solve_determinant()