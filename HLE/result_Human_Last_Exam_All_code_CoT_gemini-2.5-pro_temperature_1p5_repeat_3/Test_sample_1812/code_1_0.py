import numpy as np

def compute_determinant_and_show_steps():
    """
    Computes the determinant of a specific 3x3 matrix and shows the calculation steps.
    """
    # The given matrix from the problem
    # The context about the Markov quiver is not needed for the calculation itself.
    A = np.array([
        [2, -1, -1],
        [-1, 2, -1],
        [-1, -1, 2]
    ])

    print("Given the matrix A:")
    print(A)
    print("\nWe will compute the determinant using the formula for a 3x3 matrix:")
    print("det = a(ei - fh) - b(di - fg) + c(dh - eg)")
    
    # Extracting elements for clarity in the formula
    a, b, c = A[0, 0], A[0, 1], A[0, 2]
    d, e, f = A[1, 0], A[1, 1], A[1, 2]
    g, h, i = A[2, 0], A[2, 1], A[2, 2]

    # Printing the formula with the numbers substituted
    print("\nSubstituting the values from matrix A into the formula:")
    final_equation = (
        f"det(A) = ({a}) * (({e})*({i}) - ({f})*({h})) - "
        f"({b}) * (({d})*({i}) - ({f})*({g})) + "
        f"({c}) * (({d})*({h}) - ({e})*({g}))"
    )
    print(final_equation)

    # Calculating the intermediate parts of the equation
    part1_val = e * i - f * h
    part2_val = d * i - f * g
    part3_val = d * h - e * g
    
    step2_equation = (
        f"       = ({a}) * ({part1_val}) - ({b}) * ({part2_val}) + ({c}) * ({part3_val})"
    )
    print(step2_equation)
    
    # Calculating the final terms
    term1 = a * part1_val
    term2 = -b * part2_val
    term3 = c * part3_val

    step3_equation = f"       = {term1} + {term2} + {term3}"
    print(step3_equation)
    
    # The final result of the manual calculation
    final_result_manual = term1 + term2 + term3
    print(f"       = {final_result_manual}")
    
    # For verification, we compute using numpy's built-in function
    determinant = np.linalg.det(A)
    
    # Printing the final answer
    print("\n-------------------------------------------")
    print(f"The determinant of the matrix is: {round(determinant)}")
    print("-------------------------------------------")

# Execute the function
compute_determinant_and_show_steps()