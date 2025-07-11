import numpy as np

def calculate_brockett_minimum(a, b, s_A, s_B):
    """
    Calculates the minimum of the asymmetric Brockett cost function.

    Args:
        a (list or np.array): Singular values of matrix A, sorted descending.
        b (list or np.array): Singular values of matrix B, sorted descending.
        s_A (int): Sign of the determinant of A (-1, 0, or 1).
        s_B (int): Sign of the determinant of B (-1, 0, or 1).
    """
    n = len(a)
    if n != len(b):
        raise ValueError("Singular value lists must have the same length.")

    a = np.array(a)
    b = np.array(b)

    # b_rev corresponds to the sequence b_{n-i+1}
    b_rev = b[::-1]

    if s_A * s_B != -1:
        # Case 1: s_A * s_B is 1 or 0
        min_val = -np.sum(a * b_rev)
        
        # Print the calculation details
        print(f"Case 1: s(det(A)) * s(det(B)) = {s_A * s_B} != -1")
        print("The minimum value is -sum(a_i * b_{n-i+1}) for i=1 to n.")
        
        sum_terms_str = " + ".join([f"{a[i]}*{b_rev[i]}" for i in range(n)])
        term_values = a * b_rev
        
        print(f"min_f = -({sum_terms_str})")
        print(f"min_f = -({'+'.join(map(str, term_values))})")
        print(f"min_f = -({np.sum(term_values)})")
        print(f"min_f = {min_val}")
    else:
        # Case 2: s_A * s_B is -1
        sum_part = np.sum(a[:-1] * b_rev[:-1])
        last_term = a[-1] * b[0] # a_n * b_1
        min_val = -sum_part + last_term

        # Print the calculation details
        print(f"Case 2: s(det(A)) * s(det(B)) = {s_A * s_B} == -1")
        print("The minimum value is -sum(a_i * b_{n-i+1}) for i=1 to n-1, plus a_n * b_1.")

        sum_terms_str = " + ".join([f"{a[i]}*{b_rev[i]}" for i in range(n - 1)])
        term_values = a[:-1] * b_rev[:-1]
        
        print(f"min_f = -({sum_terms_str}) + {a[-1]}*{b[0]}")
        if len(term_values) > 0:
            print(f"min_f = -({'+'.join(map(str, term_values))}) + {last_term}")
            print(f"min_f = -({np.sum(term_values)}) + {last_term}")
        else: # Case n=1
            print(f"min_f = -(0) + {last_term}")

        print(f"min_f = {-np.sum(term_values)} + {last_term}")
        print(f"min_f = {min_val}")

# --- Example 1 ---
print("--- Example 1 ---")
# Let n=4.
a1 = [4, 3, 2, 1]
b1 = [5, 3, 2, 1]
# Assume det(A) > 0 and det(B) > 0
s_A1 = 1
s_B1 = 1
calculate_brockett_minimum(a1, b1, s_A1, s_B1)
print("\n" + "="*20 + "\n")

# --- Example 2 ---
print("--- Example 2 ---")
# Let n=4.
a2 = [4, 3, 2, 1]
b2 = [5, 3, 2, 1]
# Assume det(A) > 0 and det(B) < 0
s_A2 = 1
s_B2 = -1
calculate_brockett_minimum(a2, b2, s_A2, s_B2)
