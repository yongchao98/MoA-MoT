import collections

def continuant(coeffs):
    """
    Computes the numerator of a continued fraction (continuant) for a list of coefficients.
    K_n = a_n * K_{n-1} + K_{n-2}
    Initial conditions: K_0 = 1, K_1 = a_1
    """
    if not isinstance(coeffs, collections.abc.Sequence):
        # Handle single number input
        return coeffs
    if len(coeffs) == 0:
        return 1
    if len(coeffs) == 1:
        return coeffs[0]
    
    # Iterative calculation based on the recurrence relation
    k_minus_2 = 1
    k_minus_1 = coeffs[0]
    
    for i in range(1, len(coeffs)):
        k_current = coeffs[i] * k_minus_1 + k_minus_2
        k_minus_2 = k_minus_1
        k_minus_1 = k_current
        
    return k_minus_1

def solve_for_ck(k, a):
    """
    Solves for c_k for a given k and a list of coefficients a_1, ..., a_k.
    """
    if k < 2:
        print("Error: k must be greater than or equal to 2.")
        return

    # Ensure 'a' has at least k elements (it's 0-indexed)
    if len(a) < k:
        print(f"Error: The list 'a' must contain at least {k} elements.")
        return

    # Construct the sequence for the LHS: N[a_2,..., a_{k}+1, a_k,...,a_1]
    # Python slicing a[1:k-1] corresponds to a_2, ..., a_{k-1}
    # a[k-1] corresponds to a_k
    # a[k-1::-1] corresponds to a_k, ..., a_1
    seq_LHS = list(a[1:k-1]) + [a[k-1] + 1] + list(reversed(a[:k]))
    
    # Compute the LHS value
    val_LHS = continuant(seq_LHS)

    # Construct the sequence for the RHS term: N[a_1,...,a_{k}, a_k,...,a_2]
    # a[:k] corresponds to a_1, ..., a_k
    # a[k-1:0:-1] corresponds to a_k, ..., a_2
    seq_RHS1 = list(a[:k]) + list(a[k-1:0:-1])

    # Compute the RHS term value
    val_RHS1 = continuant(seq_RHS1)

    # Calculate c_k by subtraction
    c_k = val_LHS - val_RHS1

    # Print the final equation with the calculated numbers
    print("Verifying the equation: N[a_2,..., a_k+1, ...,a_1] = N[a_1,...,a_k,...,a_2] + c_k")
    print(f"For k={k} and a={a[:k]}:")
    print(f"{val_LHS} = {val_RHS1} + {c_k}")
    print("-" * 20)

    # Verify with the derived formula: c_k = N[a_2, ..., a_{k-1}] * N[a_1, ..., a_k]
    # N[a_2, ..., a_{k-1}] corresponds to continuant(a[1:k-1])
    # N[a_1, ..., a_k] corresponds to continuant(a[:k])
    term1_formula = continuant(a[1:k-1])
    term2_formula = continuant(a[:k])
    ck_formula = term1_formula * term2_formula
    
    print("Verifying the derived formula for c_k:")
    print(f"c_k = N[a_2,...,a_{k-1}] * N[a_1,...,a_k]")
    print(f"c_k = {term1_formula} * {term2_formula} = {ck_formula}")
    
    # Final check
    if c_k == ck_formula:
        print("The derived formula is correct for this case.")
    else:
        print("There is a mismatch in the formula.")


# --- Example Usage ---
# Define k and the coefficients a_i
k_example = 4
a_coeffs = [1, 2, 3, 4, 5] # a_1=1, a_2=2, ...

solve_for_ck(k_example, a_coeffs)
