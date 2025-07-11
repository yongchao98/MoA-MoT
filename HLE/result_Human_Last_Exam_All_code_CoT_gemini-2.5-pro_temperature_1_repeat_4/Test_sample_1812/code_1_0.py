import numpy as np

def compute_determinant_and_show_equation():
    """
    Computes the determinant of the given 3x3 matrix and prints the full calculation equation.
    """
    # The adjacency matrix is given as:
    # A = [[ 2, -1, -1],
    #      [-1,  2, -1],
    #      [-1, -1,  2]]
    A = np.array([
        [2, -1, -1],
        [-1, 2, -1],
        [-1, -1, 2]
    ])

    # Extract elements for the formula: det([[a,b,c],[d,e,f],[g,h,i]]) = a(ei-fh) - b(di-fg) + c(dh-eg)
    a, b, c = A[0, 0], A[0, 1], A[0, 2]
    d, e, f = A[1, 0], A[1, 1], A[1, 2]
    g, h, i = A[2, 0], A[2, 1], A[2, 2]

    # Calculate the determinant value
    determinant = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)

    # Construct the equation string to show the step-by-step calculation with original numbers
    # Note: We use parentheses to make the substitution of negative numbers clear.
    equation_str = (
        f"det(A) = {a} * ({e} * {i} - ({f}) * ({h})) "
        f"- ({b}) * (({d}) * {i} - ({f}) * ({g})) "
        f"+ ({c}) * (({d}) * ({h}) - {e} * ({g}))"
    )

    # Print the full equation with the final result
    print(f"The equation for the determinant is:")
    print(f"{equation_str} = {int(determinant)}")

# Execute the function
compute_determinant_and_show_equation()