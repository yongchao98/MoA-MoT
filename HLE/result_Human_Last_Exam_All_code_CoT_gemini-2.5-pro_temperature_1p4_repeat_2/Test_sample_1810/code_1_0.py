import functools

# Define a function to compute the numerator of a continued fraction K(a_1, ..., a_n)
# functools.lru_cache is used for memoization to handle recursive calls efficiently.
@functools.lru_cache(maxsize=None)
def K(a):
    """Computes the numerator of the continued fraction [a_1, ..., a_n]."""
    if len(a) == 0:
        return 1
    if len(a) == 1:
        return a[0]
    return a[-1] * K(a[:-1]) + K(a[:-2])

def solve_for_ck(k, a_list):
    """
    Solves for c_k based on the given equation.
    
    The equation is:
    N[a_2,..., a_{k}+1, a_k,...,a_1] = N[a_1,...,a_{k}, a_k,...,a_2] + c_k
    """
    if not isinstance(k, int) or k < 2:
        raise ValueError("k must be an integer greater than or equal to 2.")
    if len(a_list) != k:
        raise ValueError(f"The list 'a_list' must contain k={k} elements.")
        
    a = tuple(a_list) # Use tuple to make it hashable for the cache

    # Construct the list of coefficients for the LHS
    # seq_L = (a_2, ..., a_{k-1}, a_k+1, a_k, ..., a_1)
    # Note: Python indices are 0-based, so a_1 is a[0]
    seq_L = list(a[1:k-1]) + [a[k-1] + 1] + list(reversed(a))
    
    # Construct the list of coefficients for the RHS term
    # seq_R = (a_1, ..., a_k, a_k, ..., a_2)
    seq_R = list(a) + list(reversed(a[1:]))

    # Calculate the values
    val_lhs = K(tuple(seq_L))
    val_rhs = K(tuple(seq_R))
    
    # Calculate c_k
    ck_val = val_lhs - val_rhs

    # Print the equation with the computed values
    print(f"For k={k} and a={list(a)}:")
    print("The equation is:")
    print(f"N[{', '.join(map(str, seq_L))}] = N[{', '.join(map(str, seq_R))}] + c_k")
    print("Substituting the calculated values:")
    print(f"{val_lhs} = {val_rhs} + {ck_val}")

    # Verify with the derived formula: c_k = K(a_1..k) * K(a_2..k-1)
    p_k = K(a)
    q_k_minus_1 = K(a[1:k-1])
    ck_formula_val = p_k * q_k_minus_1
    print("\nVerification using the derived formula for c_k:")
    print(f"c_k = N[a_1,...,a_k] * N[a_2,...,a_{k-1}]")
    print(f"c_k = {p_k} * {q_k_minus_1} = {ck_formula_val}")
    if ck_val == ck_formula_val:
        print("The result is consistent.")
    else:
        print("There is an inconsistency.")

if __name__ == '__main__':
    # Example case: k=3 and a = [2, 3, 4]
    k_example = 3
    a_example = [2, 3, 4]
    solve_for_ck(k_example, a_example)
    
    print("\n" + "="*30 + "\n")

    # Example case: k=2 and a = [5, 6]
    k_example_2 = 2
    a_example_2 = [5, 6]
    solve_for_ck(k_example_2, a_example_2)