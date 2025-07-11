import functools

# Define a function to compute the numerator of a continued fraction, N(a) or K(a).
# We use memoization to handle the recursive calls efficiently.
@functools.lru_cache(maxsize=None)
def N(a):
    """
    Computes the numerator of the continued fraction [a_1, a_2, ..., a_n].
    The input 'a' is a tuple of coefficients.
    This function is also known as the continuant.
    """
    if not isinstance(a, tuple):
        a = tuple(a)
    
    n = len(a)
    if n == 0:
        return 1
    if n == 1:
        return a[0]
    
    # Recurrence relation: K(a_1,...,a_n) = a_n * K(a_1,...,a_{n-1}) + K(a_1,...,a_{n-2})
    return a[-1] * N(a[:-1]) + N(a[:-2])

def solve_for_ck(k, a_coeffs):
    """
    Solves for c_k given k and the coefficients a_1, ..., a_k.
    It verifies the derived formula for c_k by direct computation from the given equation.
    """
    if k < 2:
        raise ValueError("k must be >= 2")
    if len(a_coeffs) < k:
        raise ValueError(f"a_coeffs must contain at least k={k} elements")

    # The problem is generic for any positive integers a_i, so we use the first k coefficients provided.
    a = a_coeffs[:k]

    # Construct the sequence for the LHS: [a_2, ..., a_{k-1}, a_k+1, a_k, ..., a_1]
    # Note that Python indices are 0-based, so a_1 is a[0].
    seq_LHS = list(a[1:k-1]) + [a[k-1] + 1] + list(reversed(a))
    
    # Construct the sequence for the RHS term: [a_1, ..., a_k, a_k, ..., a_2]
    seq_RHS = list(a) + list(reversed(a[1:]))

    # Calculate the values from the equation
    val_LHS = N(tuple(seq_LHS))
    val_RHS_term = N(tuple(seq_RHS))
    
    # Solve for c_k from the equation
    c_k_from_eq = val_LHS - val_RHS_term

    print(f"For k={k} and a={a}:")
    print("The equation is:")
    print(f"N{[a_i for a_i in seq_LHS]} = N{[a_i for a_i in seq_RHS]} + c_k")
    print("\nSubstituting the calculated values:")
    # We use the instruction to output each number in the final equation.
    print(f"{val_LHS} = {val_RHS_term} + {c_k_from_eq}")

    # Verify with the derived symbolic formula: c_k = N[a_2, ..., a_{k-1}] * N[a_1, ..., a_k]
    # For k=2, the sequence a_2,...,a_{k-1} is empty, and N([])=1.
    seq_c1 = a[1:k-1]
    seq_c2 = a[:k]
    c_k_from_formula = N(tuple(seq_c1)) * N(tuple(seq_c2))
    
    print("\nVerification using the derived formula for c_k:")
    print(f"c_k = N[a_2,...,a_{k-1}] * N[a_1,...,a_k]")
    print(f"c_k = N({list(seq_c1)}) * N({list(seq_c2)})")
    print(f"c_k = {N(tuple(seq_c1))} * {N(tuple(seq_c2))} = {c_k_from_formula}")
    
    if c_k_from_eq == c_k_from_formula:
        print("\nThe results match.")
    else:
        print("\nError: The results do not match.")

# --- Example Usage ---
# Let's test for k=3 with coefficients a_1=1, a_2=2, a_3=3.
k_example = 3
a_example = [1, 2, 3, 4] # We only need the first k elements
solve_for_ck(k_example, a_example)

print("\n" + "="*30 + "\n")

# Let's test for k=2 with coefficients a_1=4, a_2=5.
k_example_2 = 2
a_example_2 = [4, 5]
solve_for_ck(k_example_2, a_example_2)
