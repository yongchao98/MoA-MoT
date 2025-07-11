import numpy as np

def solve_brockett_minimum():
    """
    Calculates the minimum of the asymmetric Brockett cost function for given matrices A and B.
    """
    # Define two 3x3 matrices A and B for demonstration.
    # You can change these matrices to any n x n matrices.
    A = np.array([
        [4, 1, 1],
        [1, 3, 2],
        [0, 1, 5]
    ])
    
    B = np.array([
        [6, 0, 1],
        [0, 2, 1],
        [1, 1, 3]
    ])
    
    n = A.shape[0]
    if A.shape != B.shape or A.shape[0] != A.shape[1]:
        print("Matrices A and B must be square and of the same size.")
        return

    # Step 1: Compute singular values of A and B.
    # np.linalg.svdvals returns them in descending order.
    a = np.linalg.svdvals(A)
    b = np.linalg.svdvals(B)
    
    print(f"Let A = \n{A}")
    print(f"Let B = \n{B}")
    print("-" * 30)
    
    # Step 2: Compute determinants and their signs.
    det_A = np.linalg.det(A)
    det_B = np.linalg.det(B)
    
    s_det_A = np.sign(det_A)
    s_det_B = np.sign(det_B)
    
    s_det_AB = np.sign(det_A * det_B)

    # Step 3: Apply the formula for the minimum value.
    sum_part = np.sum(a[:-1] * b[:-1])
    last_term = s_det_AB * a[-1] * b[-1]
    min_value = -sum_part - last_term
    
    # Step 4: Print the results and the calculation details.
    print(f"The dimension is n = {n}")
    print(f"The singular values of A are a = {a}")
    print(f"The singular values of B are b = {b}")
    print(f"The determinant of A is {det_A:.4f}, so s(det(A)) = {s_det_A:.0f}")
    print(f"The determinant of B is {det_B:.4f}, so s(det(B)) = {s_det_B:.0f}")
    print(f"s(det(A)det(B)) = {s_det_AB:.0f}")
    print("-" * 30)
    
    print("The minimum of the cost function is calculated by the formula:")
    print(f"min = -(sum_{{i=1}}^{{n-1}} a_i * b_i) - s(det(A)det(B)) * a_n * b_n")

    # Constructing the equation string
    equation_str = "-("
    for i in range(n - 1):
        equation_str += f"{a[i]:.4f} * {b[i]:.4f}"
        if i < n - 2:
            equation_str += " + "
    equation_str += f") - {s_det_AB:.0f} * {a[n-1]:.4f} * {b[n-1]:.4f} = {min_value:.4f}"
    
    print("\nCalculation:")
    print(equation_str)
    
    print(f"\nFinal calculated minimum value: {min_value:.4f}")

solve_brockett_minimum()