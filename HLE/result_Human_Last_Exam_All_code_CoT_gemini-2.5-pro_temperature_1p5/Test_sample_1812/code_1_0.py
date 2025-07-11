import numpy as np

def solve():
    """
    Computes the determinant of the given 3x3 matrix and shows the calculation steps.
    """
    # The given matrix is:
    # A = [[2, -1, -1],
    #      [-1,  2, -1],
    #      [-1, -1,  2]]

    # Assign elements to variables
    a, b, c = 2, -1, -1
    d, e, f = -1, 2, -1
    g, h, i = -1, -1, 2

    print("The given matrix A is:")
    A = np.array([[a, b, c], [d, e, f], [g, h, i]])
    print(A)
    print("\nThe determinant is calculated using the formula for a 3x3 matrix:")
    print("det(A) = a(ei - fh) - b(di - fg) + c(dh - eg)")

    print("\nStep 1: Substitute the matrix values into the formula.")
    print(f"det(A) = ({a}) * (({e})*({i}) - ({f})*({h})) - ({b}) * (({d})*({i}) - ({f})*({g})) + ({c}) * (({d})*({h}) - ({e})*({g}))")

    # Calculate the 2x2 sub-determinants
    sub_det1 = e * i - f * h
    sub_det2 = d * i - f * g
    sub_det3 = d * h - e * g

    print("\nStep 2: Calculate the values of the 2x2 determinants.")
    print(f"det(A) = ({a})*({sub_det1}) - ({b})*({sub_det2}) + ({c})*({sub_det3})")

    # Multiply by the top row elements
    term1 = a * sub_det1
    term2 = b * sub_det2
    term3 = c * sub_det3

    print("\nStep 3: Perform the multiplications.")
    print(f"det(A) = ({term1}) - ({term2}) + ({term3})")

    # Calculate the final determinant
    determinant = term1 - term2 + term3
    print("\nStep 4: Sum the terms to get the final result.")
    print(f"det(A) = {term1} - {term2} + {term3} = {determinant}")

    print("\nFinal Answer:")
    print(f"The determinant of the matrix is: {determinant}")

solve()