import math

def get_rational_cf(p, q):
    """
    Computes the continued fraction coefficients for a rational number p/q.
    Example: get_rational_cf(4, 7) -> [0, 1, 1, 3]
    """
    coeffs = []
    if q == 0:
        raise ValueError("Denominator cannot be zero")
    while q != 0:
        a = p // q
        coeffs.append(a)
        p, q = q, p % q
    return coeffs

def get_cf_matrix(coeffs):
    """
    Computes the convergents matrix for a sequence of CF coefficients.
    For a list of coefficients [c0, c1, ..., cn], this computes the matrix product
    M(c0) * M(c1) * ... * M(cn). The result is [[p_n, p_{n-1}], [q_n, q_{n-1}]].
    """
    p_k, p_k_minus_1 = 1, 0
    q_k, q_k_minus_1 = 0, 1
    
    for c in coeffs:
        p_next = c * p_k + p_k_minus_1
        q_next = c * q_k + q_k_minus_1
        
        p_k_minus_1, p_k = p_k, p_next
        q_k_minus_1, q_k = q_k, q_next
        
    return [[p_k, p_k_minus_1], [q_k, q_k_minus_1]]

def main():
    """
    Main function to compute the quadratic equation for the continued fraction
    associated with the generalized Markov number m_{4/7}.
    """
    # 1. Define the rational number
    p, q = 4, 7

    # 2. Compute its continued fraction and get the "word"
    cf_coeffs = get_rational_cf(p, q)
    word = cf_coeffs[1:]

    # 3. Form the palindromic period for the associated quadratic irrational
    period = word + list(reversed(word))

    # 4. Find the quadratic equation for the number with this periodic CF
    matrix = get_cf_matrix(period)
    p_last, p_penultimate = matrix[0][0], matrix[0][1]
    q_last, q_penultimate = matrix[1][0], matrix[1][1]

    # The quadratic equation for x = [period...] is derived from the relation
    # x = (p_last * x + p_penultimate) / (q_last * x + q_penultimate), which gives:
    # q_last * x^2 + (q_penultimate - p_last) * x - p_penultimate = 0
    A = q_last
    B = q_penultimate - p_last
    C = -p_penultimate

    # 5. Simplify the equation by dividing by GCD of coefficients
    common_divisor = math.gcd(math.gcd(A, B), C)
    A = A // common_divisor
    B = B // common_divisor
    C = C // common_divisor

    # 6. Output the result
    print(f"The continued fraction of {p}/{q} is: {cf_coeffs}")
    print(f"The generating sequence (word) is: {word}")
    print(f"The periodic part of the associated continued fraction is: {period}")
    print("\nLet x be the value of this continued fraction. It is a root of the quadratic equation Ax^2 + Bx + C = 0.")
    print("The coefficients of this equation are:")
    print(f"A = {A}")
    print(f"B = {B}")
    print(f"C = {C}")
    print("\nThe equation is:")
    print(f"{A}x^2 - {abs(B)}x - {abs(C)} = 0")
    
    print("\nThe numbers in the final equation are:")
    print(A)
    print(B)
    print(C)


if __name__ == "__main__":
    main()