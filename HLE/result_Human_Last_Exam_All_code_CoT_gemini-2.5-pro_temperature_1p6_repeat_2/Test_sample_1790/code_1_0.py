import math

def sigma_k(n, k):
    """
    Computes the sum of the k-th powers of the positive divisors of n.
    sigma_k(n) = sum_{d|n, d>0} d^k
    """
    if n == 0:
        return 0
    if n == 1:
        return 1
    
    total_sum = 0
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            total_sum += i**k
            if i*i != n:
                total_sum += (n//i)**k
    return total_sum

def get_eisenstein_coeffs(k, num_coeffs):
    """
    Computes the coefficients of the normalized Eisenstein series E_k(z).
    The formula is E_k(z) = 1 + C_k * sum_{n=1 to inf} sigma_{k-1}(n) * q^n.
    C_4 = 240, C_8 = 480.
    """
    if k == 8:
        C_k = 480
    else:
        # We only need k=8 for this problem
        raise ValueError("This implementation only supports k=8")
    
    coeffs = [0] * num_coeffs
    coeffs[0] = 1
    for n in range(1, num_coeffs):
        coeffs[n] = C_k * sigma_k(n, k - 1)
    return coeffs

def main():
    """
    Solves the problem by calculating coefficients and summing them up.
    """
    # We need the first three non-zero coefficients, which will be a_1, a_2, a_3.
    # So we need to calculate up to n=3. num_coeffs should be 4 (for n=0,1,2,3).
    num_coeffs = 4

    # 1. Get coefficients for E_8(z)
    e8_coeffs = get_eisenstein_coeffs(8, num_coeffs)

    # 2. Get coefficients for E_8(2z)
    # The q-expansion of g(dz) is sum a_n q^{dn}.
    e8_2z_coeffs = [0] * num_coeffs
    e8_2z_coeffs[0] = 1
    for n in range(1, num_coeffs):
        if n % 2 == 0:
            idx = n // 2
            if idx < len(e8_coeffs):
                # We need the coefficient a_{n/2} from E_8(z)
                 e8_2z_coeffs[n] = e8_coeffs[idx]
            else: # this case won't be hit with our num_coeffs
                 e8_2z_coeffs[n] = get_eisenstein_coeffs(8, idx+1)[idx]

    # 3. Construct the unnormalized cusp form g(z) = E_8(z) - E_8(2z)
    g_coeffs = [e8 - e8_2z for e8, e8_2z in zip(e8_coeffs, e8_2z_coeffs)]
    
    # The first coefficient (a_0) should be 0, a_1 is the first non-zero one
    # a_1 is at index 1
    normalizing_factor = g_coeffs[1]
    
    # 4. Normalize the cusp form to f(z)
    f_coeffs = [c / normalizing_factor for c in g_coeffs]
    
    # 5. The first three non-zero coefficients are at indices 1, 2, and 3
    c1 = int(f_coeffs[1])
    c2 = int(f_coeffs[2])
    c3 = int(f_coeffs[3])
    
    total_sum = c1 + c2 + c3
    
    # 6. Print the result as requested
    print(f"The unique normalized cusp form is f(z) = {c1}q + {c2}q^2 + {c3}q^3 + ...")
    print(f"The sum of the first three non-zero coefficients is:")
    print(f"{c1} + {c2} + {c3} = {total_sum}")

if __name__ == "__main__":
    main()
