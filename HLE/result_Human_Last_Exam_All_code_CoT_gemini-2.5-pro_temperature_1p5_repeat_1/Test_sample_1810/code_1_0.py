def get_numerator(coeffs):
    """
    Calculates the numerator of a continued fraction [a_1, a_2, ..., a_n].
    The numerator p_n is calculated using the recurrence relation:
    p_n = a_n * p_{n-1} + p_{n-2}, with initial conditions p_0 = 1 and p_{-1} = 0.
    """
    p_minus_2 = 0
    p_minus_1 = 1
    for x in coeffs:
        p = x * p_minus_1 + p_minus_2
        p_minus_2 = p_minus_1
        p_minus_1 = p
    return p

def solve_for_ck_case(k, a):
    """
    Solves for c_k for a specific case of k and coefficients a_i.
    a should be a list of k integers for a_1, ..., a_k.
    Note: Python lists are 0-indexed, so a_i corresponds to a[i-1].
    """
    if k < 2:
        raise ValueError("k must be greater than or equal to 2.")
    if len(a) != k:
        raise ValueError(f"The list of coefficients 'a' must have {k} elements.")

    # Construct the coefficient list for the LHS: N[a_2,..., a_{k}+1, a_k,...,a_1]
    lhs_coeffs = []
    # a_2, ..., a_k
    for i in range(1, k): # from a_2 to a_k
        lhs_coeffs.append(a[i])
    # a_k + 1
    lhs_coeffs.append(a[k-1] + 1)
    # a_k, ..., a_1
    for i in range(k - 1, -1, -1): # from a_k to a_1
        lhs_coeffs.append(a[i])

    # Construct the coefficient list for the RHS: N[a_1,...,a_{k}, a_k,...,a_2]
    rhs_coeffs = []
    # a_1, ..., a_k
    for i in range(k):
        rhs_coeffs.append(a[i])
    # a_k, ..., a_2
    for i in range(k - 1, 0, -1): # from a_k to a_2
        rhs_coeffs.append(a[i])

    lhs_numerator = get_numerator(lhs_coeffs)
    rhs_numerator = get_numerator(rhs_coeffs)
    
    c_k = lhs_numerator - rhs_numerator
    
    print(f"For k={k} and a={a}:")
    print(f"LHS Numerator N({lhs_coeffs}) = {lhs_numerator}")
    print(f"RHS Numerator N({rhs_coeffs}) = {rhs_numerator}")
    print("The equation is:")
    print(f"{lhs_numerator} = {rhs_numerator} + c_{k}")
    print(f"Solving for c_{k}:")
    print(f"c_{k} = {lhs_numerator} - {rhs_numerator} = {c_k}")


# Example case: k=2, with coefficients a_1=1, a_2=2
k_val = 2
a_val = [1, 2]
solve_for_ck_case(k_val, a_val)