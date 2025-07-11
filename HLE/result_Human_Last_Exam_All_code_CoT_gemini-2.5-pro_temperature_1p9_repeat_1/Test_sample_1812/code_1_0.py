import numpy as np

def compute_determinant():
    """
    Computes the determinant of the given 3x3 matrix and shows the calculation steps.
    """
    # The given adjacency matrix A
    A = np.array([
        [2, -1, -1],
        [-1, 2, -1],
        [-1, -1, 2]
    ])

    print("The given matrix is:")
    print(A)
    print("\nTo compute the determinant, we use the formula for a 3x3 matrix:")
    print("det(A) = a(ei - fh) - b(di - fg) + c(dh - eg)")

    # Extract elements for clarity in the equation
    a, b, c = A[0, 0], A[0, 1], A[0, 2]
    d, e, f = A[1, 0], A[1, 1], A[1, 2]
    g, h, i = A[2, 0], A[2, 1], A[2, 2]

    # Step 1: Substitute the numbers from the matrix into the formula
    equation_step1 = (f"det(A) = ({a}) * (({e})*({i}) - ({f})*({h}))"
                      f" - ({b}) * (({d})*({i}) - ({f})*({g}))"
                      f" + ({c}) * (({d})*({h}) - ({e})*({g}))")
    print("\nStep 1: Substituting the numbers into the formula gives:")
    print(equation_step1)

    # Step 2: Simplify the expressions inside the parentheses
    minor1 = e*i - f*h
    minor2 = d*i - f*g
    minor3 = d*h - e*g
    
    equation_step2 = f"det(A) = ({a})*({minor1}) - ({b})*({minor2}) + ({c})*({minor3})"
    print("\nStep 2: Simplifying the 2x2 determinants inside the parentheses:")
    print(equation_step2)

    # Step 3: Perform the multiplications
    term1 = a * minor1
    term2 = -b * minor2
    term3 = c * minor3
    
    equation_step3 = f"det(A) = {term1} + ({term2}) + ({term3})"
    print("\nStep 3: Performing the multiplications:")
    print(equation_step3)
    
    # Step 4: Calculate the final determinant
    determinant = np.linalg.det(A)
    
    # We use round() to handle potential floating point inaccuracies for a result that should be an integer
    final_determinant = round(determinant)
    
    print(f"\nStep 4: Summing the terms gives the final determinant:")
    print(f"det(A) = {final_determinant}")

compute_determinant()
<<<0>>>