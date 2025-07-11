import numpy as np

def solve():
    """
    Computes the determinant of the given adjacency matrix A.
    """
    # Define the matrix A
    A = np.array([
        [2, -1, -1],
        [-1, 2, -1],
        [-1, -1, 2]
    ])

    # Extract elements for the formula
    a, b, c = A[0]
    d, e, f = A[1]
    g, h, i = A[2]

    # Calculate the determinant using the formula
    # det = a(ei - fh) - b(di - fg) + c(dh - eg)
    term1 = a * (e * i - f * h)
    term2 = -b * (d * i - f * g)
    term3 = c * (d * h - e * g)
    determinant = term1 + term2 + term3

    print("The given adjacency matrix is:")
    print(A)
    print("\nThe formula for the determinant is: det(A) = a(ei - fh) - b(di - fg) + c(dh - eg)")
    print("\nSubstituting the values from the matrix:")
    
    # Print the equation with numbers
    # We construct the string carefully to show the equation clearly
    # det = 2 * (2*2 - (-1)*(-1)) - (-1) * ((-1)*2 - (-1)*(-1)) + (-1) * ((-1)*(-1) - 2*(-1))
    print("det(A) = {}*({}*{} - {}*{}) - ({})*({}*{} - {}*{}) + ({})*({}*{} - {}*{})".format(
        a, e, i, f, h, b, d, i, f, g, c, d, h, e, g
    ))

    # Print intermediate calculation
    # det = 2 * (4 - 1) + 1 * (-2 - 1) - 1 * (1 - -2)
    print("det(A) = {}*({}) - ({})*({}) + ({})*({})".format(
        a, e*i - f*h, b, d*i - f*g, c, d*h - e*g
    ))
    
    # Print the final steps of the arithmetic
    # det = 2*3 + 1*(-3) - 1*3
    print("det(A) = {} + ({}) + ({})".format(term1, term2, term3))
    print(f"det(A) = {determinant}")

solve()
<<<0>>>