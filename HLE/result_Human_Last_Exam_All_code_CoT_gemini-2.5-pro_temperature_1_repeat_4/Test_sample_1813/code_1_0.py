import math

def get_cf_coeffs(p, q):
    """
    Computes the coefficients of the continued fraction of the rational number p/q.
    """
    coeffs = []
    if q == 0:
        return []
    while q != 0:
        a = p // q
        coeffs.append(a)
        p, q = q, p % q
    return coeffs

def multiply_matrices(A, B):
    """
    Multiplies two 2x2 matrices A and B.
    """
    C = [[0, 0], [0, 0]]
    C[0][0] = A[0][0] * B[0][0] + A[0][1] * B[1][0]
    C[0][1] = A[0][0] * B[0][1] + A[0][1] * B[1][1]
    C[1][0] = A[1][0] * B[0][0] + A[1][1] * B[1][0]
    C[1][1] = A[1][0] * B[0][1] + A[1][1] * B[1][1]
    return C

def main():
    """
    Main function to compute the continued fraction for the generalized Markov number m_{4/7}.
    """
    # The rational number is p/q = 4/7.
    p, q = 4, 7

    # Step 1: Compute the continued fraction of 4/7.
    # For a rational r in (0,1), the generalized Markov number is defined
    # based on the continued fraction [a_1, a_2, ..., a_n], where r = [0; a_1, ...].
    cf_coeffs_pq = get_cf_coeffs(p, q)
    # We drop the initial a_0 = 0 from [0; 1, 1, 3].
    markov_cf_coeffs = cf_coeffs_pq[1:]

    # Step 2: Define and compute the matrix product M = M_1 * M_2 * ... * M_n
    # where M_i = [[a_i, 1], [1, 0]].
    if not markov_cf_coeffs:
        print("No coefficients found for matrix product.")
        return
        
    M_prod = [[markov_cf_coeffs[0], 1], [1, 0]]
    for i in range(1, len(markov_cf_coeffs)):
        a = markov_cf_coeffs[i]
        M_i = [[a, 1], [1, 0]]
        M_prod = multiply_matrices(M_prod, M_i)

    # Step 3: Compute the trace of the resulting matrix.
    trace = M_prod[0][0] + M_prod[1][1]

    # Step 4: Compute the generalized Markov number m_{p/q} = trace / 3.
    # We'll keep it as a fraction and reduce it.
    m_num = trace
    m_den = 3
    common_divisor = math.gcd(m_num, m_den)
    m_num //= common_divisor
    m_den //= common_divisor

    # Step 5: Compute the continued fraction of the resulting m_{p/q}.
    m_cf_coeffs = get_cf_coeffs(m_num, m_den)

    # Step 6: Format the output to show the final equation with all its numbers.
    if len(m_cf_coeffs) == 1:
        cf_string = f"[{m_cf_coeffs[0]}]"
    else:
        # Standard notation [a_0; a_1, a_2, ...]
        cf_string = f"[{m_cf_coeffs[0]}; {', '.join(map(str, m_cf_coeffs[1:]))}]"

    print(f"The generalized Markov number m_{p}/{q} is {m_num}/{m_den}.")
    print(f"The continued fraction of {m_num}/{m_den} is: {cf_string}")

if __name__ == "__main__":
    main()