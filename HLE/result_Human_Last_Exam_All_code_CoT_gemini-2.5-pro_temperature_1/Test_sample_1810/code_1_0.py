import random

def N(a):
    """
    Computes the numerator of the continued fraction [a_1, a_2, ...].
    a: a list of integers.
    """
    if not a:
        # By convention, the numerator of an empty continued fraction is 1.
        return 1
    if len(a) == 1:
        return a[0]
    
    # Using the recurrence p_n = a_n * p_{n-1} + p_{n-2}
    # with initial values p_0 = 1, p_1 = a_1.
    # We can compute iteratively.
    p_prev = 1
    p_curr = a[0]
    for i in range(1, len(a)):
        p_next = a[i] * p_curr + p_prev
        p_prev = p_curr
        p_curr = p_next
    return p_curr

def solve_and_verify(k, a_coeffs):
    """
    Solves for c_k and verifies the formula for a given k and coefficients.
    k: an integer >= 2
    a_coeffs: a list of k positive integers [a_1, ..., a_k]
    """
    if k < 2 or len(a_coeffs) != k:
        print("Error: k must be >= 2 and the length of a_coeffs must be k.")
        return

    print(f"--- Verifying for k={k} and a={a_coeffs} ---")
    
    # Construct the sequences for LHS and RHS
    # S_L = [a_2, ..., a_{k-1}, a_k+1, a_k, ..., a_1]
    S_L = a_coeffs[1:k-1] + [a_coeffs[k-1] + 1] + list(reversed(a_coeffs))
    
    # S_R = [a_1, ..., a_k, a_k, ..., a_2]
    S_R = a_coeffs + list(reversed(a_coeffs[1:]))

    # Calculate the values
    lhs_val = N(S_L)
    rhs_n_val = N(S_R)
    
    # Calculate c_k from the original equation
    c_k_from_eq = lhs_val - rhs_n_val
    
    # The equation is: N(S_L) = N(S_R) + c_k
    print("Original equation values:")
    print(f"N({S_L}) = {lhs_val}")
    print(f"N({S_R}) = {rhs_n_val}")
    print(f"Calculated c_{k} = {lhs_val} - {rhs_n_val} = {c_k_from_eq}")
    print(f"Final equation with numbers: {lhs_val} = {rhs_n_val} + {c_k_from_eq}")
    
    # Verify the derived formula for c_k
    # c_k = N[a_1, ..., a_k] * N[a_2, ..., a_{k-1}]
    term1 = N(a_coeffs)
    term2_seq = a_coeffs[1:k-1]
    term2 = N(term2_seq)
    c_k_from_formula = term1 * term2
    
    print("\nVerifying the derived formula for c_k:")
    print(f"c_k = N[a_1,...,a_k] * N[a_2,...,a_{k-1}]")
    print(f"N({a_coeffs}) = {term1}")
    print(f"N({term2_seq}) = {term2}")
    print(f"c_{k} from formula = {term1} * {term2} = {c_k_from_formula}")

    if c_k_from_eq == c_k_from_formula:
        print("\nSuccess! The derived formula for c_k is correct.")
    else:
        print("\nFailure! There is a discrepancy in the formula for c_k.")
        
    print("-" * (10 + len(str(k)) + len(str(a_coeffs))))

# --- Main execution ---
# You can change k and the coefficients to test other cases.
# Case 1: k=3, a=[2, 3, 4]
k_test = 3
a_test = [2, 3, 4]
solve_and_verify(k_test, a_test)

# Case 2: k=2, a=[5, 6] (to test the edge case)
k_test = 2
a_test = [5, 6]
solve_and_verify(k_test, a_test)

# Case 3: k=4, random coefficients
k_test = 4
a_test = [random.randint(1, 5) for _ in range(k_test)]
solve_and_verify(k_test, a_test)

print("\nThe general solution for c_k is:")
print("c_k = N[a_1, ..., a_k] * N[a_2, ..., a_{k-1}]")
