import collections

# Using a cache for memoization to speed up the recursive calculation for N
cache = {}

def N(a):
    """
    Computes the numerator of the continued fraction [a_1, a_2, ..., a_k].
    This is also known as the continuant K(a_1, ..., a_k).
    The input 'a' is a list of integers.
    """
    a = tuple(a)
    if a in cache:
        return cache[a]
    
    if not a:
        return 1
    if len(a) == 1:
        return a[0]
    
    # Using the recurrence relation K(a_1,...,a_k) = a_k * K(a_1,...,a_{k-1}) + K(a_1,...,a_{k-2})
    res = a[-1] * N(a[:-1]) + N(a[:-2])
    cache[a] = res
    return res

def solve_for_ck(a, k):
    """
    Solves for c_k given a list of integers 'a' and the integer 'k'.
    'a' is expected to be [a_1, a_2, ..., a_k].
    """
    if not (isinstance(a, list) and len(a) == k and k >= 2):
        print("Invalid input: 'a' must be a list of length k, and k must be >= 2.")
        return

    # Construct the sequences for the LHS and RHS of the given equation
    # Note: Python uses 0-based indexing, so a_i corresponds to a[i-1]
    
    # Sequence for LHS: [a_2, ..., a_{k-1}, a_k+1, a_k, ..., a_1]
    a_LHS = a[1:k-1] + [a[k-1] + 1] + list(reversed(a))
    
    # Sequence for RHS: [a_1, ..., a_k, a_k, ..., a_2]
    a_RHS = a + list(reversed(a[1:]))
    
    # Calculate the numerators
    N_LHS = N(a_LHS)
    N_RHS = N(a_RHS)
    
    # Calculate c_k from the equation
    c_k_from_eq = N_LHS - N_RHS
    
    # The equation is: N[LHS] = N[RHS] + c_k
    print(f"For k={k} and a={a}:")
    print(f"The given equation is:")
    print(f"N{[a_2, ..., a_k+1, a_k, ..., a_1]} = N{[a_1, ..., a_k, a_k, ..., a_2]} + c_{k}")
    print(f"Numerically:")
    print(f"{N_LHS} = {N_RHS} + c_{k}")
    print(f"So, c_{k} = {N_LHS} - {N_RHS} = {c_k_from_eq}\n")
    
    # Now, we verify this with our derived formula:
    # c_k = N[a_1, ..., a_k] * N[a_2, ..., a_{k-1}]
    
    # Sequence for the first term of the formula
    a_formula_1 = a
    # Sequence for the second term of the formula
    a_formula_2 = a[1:k-1]
    
    # Calculate the terms
    N_formula_1 = N(a_formula_1)
    N_formula_2 = N(a_formula_2)
    
    # Calculate c_k from the formula
    c_k_from_formula = N_formula_1 * N_formula_2
    
    print(f"Our derived formula is:")
    print(f"c_k = N[a_1, ..., a_k] * N[a_2, ..., a_{k-1}]")
    print(f"Numerically:")
    print(f"N[{a_formula_1}] = {N_formula_1}")
    print(f"N[{a_formula_2}] = {N_formula_2}")
    print(f"c_{k} = {N_formula_1} * {N_formula_2} = {c_k_from_formula}\n")

    # Final check
    if c_k_from_eq == c_k_from_formula:
        print("The derived formula is correct for this case.")
    else:
        print("There is a mismatch. Please check the derivation or code.")


# Example case: k=3, a = [4, 3, 2]
k_example = 3
a_example = [4, 3, 2]
solve_for_ck(a_example, k_example)

print("-" * 20)

# Example case: k=2, a = [5, 2]
k_example_2 = 2
a_example_2 = [5, 2]
solve_for_ck(a_example_2, k_example_2)
