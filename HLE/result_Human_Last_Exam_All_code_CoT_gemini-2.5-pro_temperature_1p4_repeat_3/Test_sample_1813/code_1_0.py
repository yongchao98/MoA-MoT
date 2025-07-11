import math

def get_cf_coeffs(p, q):
    """Computes the continued fraction coefficients for a rational number p/q."""
    coeffs = []
    if q == 0:
        return []
    while q != 0:
        coeffs.append(p // q)
        p, q = q, p % q
    return coeffs

def matrix_multiply(A, B):
    """Multiplies two 2x2 matrices."""
    C = [[0, 0], [0, 0]]
    C[0][0] = A[0][0] * B[0][0] + A[0][1] * B[1][0]
    C[0][1] = A[0][0] * B[0][1] + A[0][1] * B[1][1]
    C[1][0] = A[1][0] * B[0][0] + A[1][1] * B[1][0]
    C[1][1] = A[1][0] * B[0][1] + A[1][1] * B[1][1]
    return C

def format_cf_equation(coeffs):
    """Formats a list of continued fraction coefficients into a string equation."""
    if not coeffs:
        return "0"
    if len(coeffs) == 1:
        return str(coeffs[0])
    
    # Recursively build the string from the tail
    def build_string(c_list):
        if len(c_list) == 1:
            return str(c_list[0])
        else:
            return f"{c_list[0]} + 1/({build_string(c_list[1:])})"

    return build_string(coeffs)

def compute_generalized_markov_cf(p, q):
    """
    Computes the continued fraction associated with the generalized Markov number m_{p/q}.
    """
    # Step 1: Get the continued fraction of p/q = [0; c1, c2, ...]
    # We only need the coefficients after the first one (the integer part, which is 0).
    cf_pq_coeffs = get_cf_coeffs(p, q)[1:]

    # Step 2: Compute the matrix product K(C)
    # K(C) = M_c1 * M_c2 * ... * M_ck
    K_C = [[1, 0], [0, 1]]  # Start with identity matrix
    for c in cf_pq_coeffs:
        M_c = [[c, 1], [1, 0]]
        K_C = matrix_multiply(K_C, M_c)

    # Step 3: Compute the generalized Markov number m_{p/q} = tr(K(C)) / q
    trace = K_C[0][0] + K_C[1][1]
    m_num = trace
    m_den = q
    
    # Simplify the fraction
    common_divisor = math.gcd(m_num, m_den)
    m_num //= common_divisor
    m_den //= common_divisor

    # Step 4: Compute the continued fraction of m_{p/q}
    m_cf_coeffs = get_cf_coeffs(m_num, m_den)

    # Step 5: Format the output
    cf_eq_str = format_cf_equation(m_cf_coeffs)

    print(f"The continued fraction for p/q = {p}/{q} is [0; {', '.join(map(str, cf_pq_coeffs))}].")
    print(f"The matrix K(C) is [[{K_C[0][0]}, {K_C[0][1]}], [{K_C[1][0]}, {K_C[1][1]}]].")
    print(f"The trace of K(C) is {trace}.")
    print(f"The generalized Markov number m_{p/q} is trace/q = {trace}/{q}.")
    print("\nResult:")
    print(f"m_{p}/{q} = {m_num}/{m_den}")
    # The prompt requests the final equation with all numbers.
    final_eq = f"{m_num}/{m_den} = {cf_eq_str}"
    
    # Let's extract and print the numbers in the equation for clarity
    all_numbers = [m_num, m_den] + m_cf_coeffs
    print(f"The final equation is: {final_eq}")
    print(f"The numbers in the final equation are: {', '.join(map(str, all_numbers))}")


# Given p/q = 4/7
p_val = 4
q_val = 7
compute_generalized_markov_cf(p_val, q_val)
