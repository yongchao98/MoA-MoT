import numpy as np

def compute_determinant():
    """
    Computes and prints the step-by-step determinant calculation 
    for the given 3x3 matrix.
    """
    # The matrix A is the adjacency matrix of the Markov quiver.
    A = np.array([
        [2, -1, -1],
        [-1, 2, -1],
        [-1, -1, 2]
    ])

    # Extract elements to display them in the formula for clarity.
    a, b, c = A[0, 0], A[0, 1], A[0, 2]
    d, e, f = A[1, 0], A[1, 1], A[1, 2]
    g, h, i = A[2, 0], A[2, 1], A[2, 2]

    print(f"The given matrix A is:\n{A}\n")

    # The user wants to see the numbers in the final equation. We will build it step-by-step.
    print("The determinant is calculated using the standard formula for a 3x3 matrix:")
    print("det(A) = a(ei - fh) - b(di - fg) + c(dh - eg)\n")

    # Step 1: Substitute the numbers into the formula.
    print("Step 1: Substitute the values from matrix A into the formula.")
    print(f"det(A) = ({a}) * (({e})*({i}) - ({f})*({h})) - ({b}) * (({d})*({i}) - ({f})*({g})) + ({c}) * (({d})*({h}) - ({e})*({g}))\n")

    # Step 2: Evaluate the expressions inside the parentheses (the minors).
    minor1 = e * i - f * h
    minor2 = d * i - f * g
    minor3 = d * h - e * g
    print("Step 2: Evaluate the expressions inside the parentheses.")
    print(f"det(A) = ({a})*({minor1}) - ({b})*({minor2}) + ({c})*({minor3})\n")
    
    # Step 3: Perform the multiplications.
    term1 = a * minor1
    # Note the formula has -b, so we apply the sign here.
    term2 = -b * minor2
    term3 = c * minor3
    print("Step 3: Perform the multiplications to get the terms of the sum.")
    print(f"det(A) = {term1} + ({term2}) + ({term3})\n")
    
    # Step 4: Sum the terms to find the final determinant.
    determinant = term1 + term2 + term3
    print("Step 4: Sum the terms to find the final result.")
    print(f"det(A) = {determinant}\n")

    print(f"Thus, the determinant of the adjacency matrix is {determinant}.")

if __name__ == '__main__':
    compute_determinant()