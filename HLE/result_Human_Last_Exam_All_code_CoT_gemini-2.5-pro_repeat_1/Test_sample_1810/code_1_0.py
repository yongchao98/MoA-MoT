def calculate_numerator(coeffs):
    """
    Calculates the numerator of a continued fraction given its coefficients.
    This implements the recurrence for continuant polynomials K(a_1, ..., a_m).
    N[c_1, ..., c_m] = K(c_1, ..., c_m)
    """
    if not coeffs:
        return 1  # K() = 1
    if len(coeffs) == 1:
        return coeffs[0] # K(a_1) = a_1

    # Recurrence: p_i = a_i * p_{i-1} + p_{i-2}
    # We use p_prev for p_{i-2} and p_curr for p_{i-1}
    p_prev = 1  # Represents K()
    p_curr = coeffs[0] # Represents K(a_1)

    for i in range(1, len(coeffs)):
        p_next = coeffs[i] * p_curr + p_prev
        p_prev = p_curr
        p_curr = p_next
    return p_curr

def solve_and_verify(k, a):
    """
    Solves for c_k and verifies the original equation for a given k and list a.
    """
    if k < 2 or len(a) < k:
        print("Error: k must be >= 2 and len(a) must be >= k.")
        return

    print(f"Verifying for k={k} and a={a[:k]}")

    # 1. Calculate the LHS numerator: N[a_2,..., a_k+1, a_k,...,a_1]
    # Sequence is (a_2, ..., a_{k-1}, a_k+1) + (a_k, ..., a_1)
    # Note: Python uses 0-based indexing for list 'a'.
    # a_1 is a[0], a_k is a[k-1]
    seq_lhs_part1 = a[1:k-1] + [a[k-1] + 1]
    seq_lhs_part2 = list(reversed(a[0:k]))
    seq_lhs = seq_lhs_part1 + seq_lhs_part2
    val_lhs = calculate_numerator(seq_lhs)

    # 2. Calculate the RHS numerator: N[a_1,...,a_k, a_k,...,a_2]
    # Sequence is (a_1, ..., a_k) + (a_k, ..., a_2)
    seq_rhs_part1 = a[0:k]
    seq_rhs_part2 = list(reversed(a[1:k]))
    seq_rhs = seq_rhs_part1 + seq_rhs_part2
    val_rhs = calculate_numerator(seq_rhs)

    # 3. Calculate c_k using the derived formula:
    # c_k = N[a_1, ..., a_k] * N[a_2, ..., a_{k-1}]
    ck_term1_seq = a[0:k]
    ck_term2_seq = a[1:k-1]
    ck_term1 = calculate_numerator(ck_term1_seq)
    ck_term2 = calculate_numerator(ck_term2_seq)
    val_ck = ck_term1 * ck_term2

    # 4. Print the final verified equation with the calculated numbers
    print("\nOriginal Equation:")
    print("N[a_2,..., a_k+1, a_k,...,a_1] = N[a_1,...,a_k, a_k,...,a_2] + c_k")
    
    print("\nVerification with numbers:")
    print(f"{val_lhs} = {val_rhs} + {val_ck}")

    if val_lhs == val_rhs + val_ck:
        print("The equation holds true.")
    else:
        print("The equation does not hold. There might be an error in the derivation.")


# --- Example Usage ---
# Let k=4 and a_i be positive integers.
k_example = 4
a_coeffs_example = [1, 2, 3, 4, 5] # We only need the first k elements

solve_and_verify(k_example, a_coeffs_example)
