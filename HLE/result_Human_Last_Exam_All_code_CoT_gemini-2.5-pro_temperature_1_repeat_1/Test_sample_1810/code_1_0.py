def N(coeffs):
    """
    Calculates the numerator of a continued fraction for a given list of coefficients.
    N([]) = 1
    N([a_1]) = a_1
    N([a_1, a_2]) = a_1*a_2 + 1
    In general, p_n = a_n * p_{n-1} + p_{n-2}
    with p_0 = 1, p_{-1} = 0.
    """
    if not coeffs:
        return 1
    elif len(coeffs) == 1:
        return coeffs[0]
    
    p_prev2 = 1  # Corresponds to p_{n-2}, starts with p_0
    p_prev1 = coeffs[0]  # Corresponds to p_{n-1}, starts with p_1
    
    for i in range(1, len(coeffs)):
        p_curr = coeffs[i] * p_prev1 + p_prev2
        p_prev2 = p_prev1
        p_prev1 = p_curr
        
    return p_prev1

def solve_c_k(a):
    """
    Solves for c_k for a given list of coefficients a = [a_1, ..., a_k].
    c_k = N([a_1, ..., a_k]) * N([a_2, ..., a_{k-1}])
    """
    k = len(a)
    if k < 2:
        print("Error: k must be >= 2")
        return

    # Calculate N([a_1, ..., a_k])
    n1 = N(a)
    
    # Calculate N([a_2, ..., a_{k-1}])
    # For k=2, this is an empty list.
    coeffs_n2 = a[1:k-1]
    n2 = N(coeffs_n2)
    
    # Calculate c_k
    c_k = n1 * n2
    
    print(f"For k = {k} and a = {a}:")
    print(f"The first term is N({a}) = {n1}")
    print(f"The second term is N({coeffs_n2}) = {n2}")
    print(f"c_{k} = N({a}) * N({coeffs_n2})")
    print(f"c_{k} = {n1} * {n2} = {c_k}")

# Example usage with k=4 and a = [2, 3, 4, 5]
# You can change these values to test other cases.
a_coeffs = [2, 3, 4, 5]
solve_c_k(a_coeffs)

# Example for k=2
print("\n" + "="*20 + "\n")
solve_c_k([6, 7])