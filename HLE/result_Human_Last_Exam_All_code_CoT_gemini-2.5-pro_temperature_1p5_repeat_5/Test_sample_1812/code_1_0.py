import numpy as np

def solve_determinant():
    """
    Computes the determinant of the given 3x3 matrix and shows the calculation steps.
    """
    # The given adjacency matrix
    A = np.array([
        [2, -1, -1],
        [-1, 2, -1],
        [-1, -1, 2]
    ])

    # Extract elements for clarity in the formula display
    a, b, c = A[0, :]
    d, e, f = A[1, :]
    g, h, i = A[2, :]

    # Calculate intermediate values for the cofactors
    cofactor_a = e * i - f * h
    cofactor_b = d * i - f * g
    cofactor_c = d * h - e * g

    # Calculate the full determinant
    determinant = a * cofactor_a - b * cofactor_b + c * cofactor_c

    print("The given adjacency matrix is:")
    print(A)
    print("\nThe formula for the determinant of a 3x3 matrix is:")
    print("det(A) = a(ei - fh) - b(di - fg) + c(dh - eg)\n")
    print("Step 1: Substitute the values from matrix A into the formula.")
    # Helper to format negative numbers with parentheses
    def s(val):
        return f"({val})" if val < 0 else str(val)

    print(f"det(A) = {s(a)}*({s(e)}*{s(i)} - {s(f)}*{s(h)}) - {s(b)}*({s(d)}*{s(i)} - {s(f)}*{s(g)}) + {s(c)}*({s(d)}*{s(h)} - {s(e)}*{s(g)})")

    print("\nStep 2: Calculate the products inside the parentheses.")
    print(f"det(A) = {s(a)}*({e*i} - {f*h}) - {s(b)}*({d*i} - {f*g}) + {s(c)}*({d*h} - {e*g})")

    print("\nStep 3: Calculate the value of each sub-determinant (cofactor term).")
    print(f"det(A) = {s(a)}*({cofactor_a}) - {s(b)}*({cofactor_b}) + {s(c)}*({cofactor_c})")

    print("\nStep 4: Perform the multiplications.")
    print(f"det(A) = {a * cofactor_a} - ({b * cofactor_b}) + ({c * cofactor_c})")

    print("\nStep 5: Perform the final addition and subtraction to get the result.")
    print(f"det(A) = {determinant}")

    print(f"\nTherefore, the determinant of the matrix A is {int(round(determinant))}.")

solve_determinant()