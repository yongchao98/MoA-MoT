import sys

def numerator(coeffs):
    """
    Computes the numerator of a finite continued fraction [x_1, ..., x_m].
    The recurrence relation is p_m = x_m * p_{m-1} + p_{m-2}.
    p_0 = 1, p_1 = x_1.
    """
    if not coeffs:
        # Numerator of an empty continued fraction is defined as 1.
        # This corresponds to p_0 in the recurrence p_2 = x_2*p_1 + p_0.
        return 1
    if len(coeffs) == 1:
        return coeffs[0]

    p_prev2 = 1
    p_prev1 = coeffs[0]

    for i in range(1, len(coeffs)):
        p_curr = coeffs[i] * p_prev1 + p_prev2
        p_prev2 = p_prev1
        p_prev1 = p_curr
    return p_prev1

def solve_for_ck(k, a):
    """
    Solves for c_k for a given k and a list of coefficients a = [a_1, ..., a_k].
    """
    if not isinstance(a, list) or len(a) != k:
        print(f"Error: The list 'a' must contain k={k} elements.", file=sys.stderr)
        return
    if k < 2:
        print(f"Error: k must be >= 2.", file=sys.stderr)
        return

    # Construct the sequence for the LHS: [a_2, ..., a_{k-1}, a_k+1, a_k, ..., a_1]
    # Note: Python indices are 0-based, so a_1 is a[0].
    seq_lhs = a[1:k-1] + [a[k-1] + 1] + a[::-1][0:k]
    
    # Construct the sequence for the RHS: [a_1, ..., a_k, a_k, ..., a_2]
    seq_rhs = a[0:k] + a[1:k][::-1]

    # Compute the numerators
    n_lhs = numerator(seq_lhs)
    n_rhs1 = numerator(seq_rhs)

    # Calculate c_k
    c_k = n_lhs - n_rhs1

    # Print the equation with the calculated values
    print(f"For k={k} and a={a}:")
    print(f"The equation is: N[{', '.join(map(str, seq_lhs))}] = N[{', '.join(map(str, seq_rhs))}] + c_k")
    print(f"Numerically: {n_lhs} = {n_rhs1} + {c_k}")
    
    # Verification using the derived formula: c_k = N[a_2, ..., a_{k-1}] * N[a_1, ..., a_k]
    # For k=2, the first term N[a_2, ..., a_{k-1}] is N[] which is 1.
    n_c1 = numerator(a[1:k-1]) 
    n_c2 = numerator(a)
    c_k_formula = n_c1 * n_c2
    print(f"Verification using formula c_k = N[a_2..a_{k-1}] * N[a_1..a_k]: {c_k_formula}")
    print("-" * 20)


# Example usage:
# Case 1: k=2, a = [1, 2]
solve_for_ck(2, [1, 2])

# Case 2: k=3, a = [1, 1, 1] (Fibonacci case)
solve_for_ck(3, [1, 1, 1])

# Case 3: k=4, a = [1, 2, 1, 3]
solve_for_ck(4, [1, 2, 1, 3])
