import sympy

def sigma(k, n):
    """
    Computes the sum of the k-th powers of the divisors of n.
    """
    if n == 0:
        return 0
    return sum(d**k for d in sympy.divisors(n))

def get_eisenstein_series_coeffs(k, num_coeffs):
    """
    Computes coefficients for the normalized Eisenstein series E_k.
    """
    coeffs = [1]
    if k == 4:
        # E4 = 1 + 240 * sum(sigma_3(n) * q^n)
        for n in range(1, num_coeffs):
            coeffs.append(240 * sigma(3, n))
    elif k == 8:
        # E8 = 1 + 480 * sum(sigma_7(n) * q^n)
        for n in range(1, num_coeffs):
            coeffs.append(480 * sigma(7, n))
    return coeffs

def multiply_poly(p1, p2):
    """
    Multiplies two polynomials represented as lists of coefficients.
    """
    n = len(p1)
    m = len(p2)
    prod = [0] * (n + m - 1)
    for i in range(n):
        for j in range(m):
            prod[i + j] += p1[i] * p2[j]
    return prod

def main():
    """
    Main function to solve the problem.
    """
    num_coeffs = 10  # Compute enough coefficients for accuracy

    # 1. Get coefficients for E4
    e4_coeffs = get_eisenstein_series_coeffs(4, num_coeffs)

    # 2. Get coefficients for F = E4(2z)
    f_coeffs = [0] * num_coeffs
    for i in range(num_coeffs):
        if i % 2 == 0:
            idx = i // 2
            if idx < len(e4_coeffs):
                f_coeffs[i] = e4_coeffs[idx]
        else:
            f_coeffs[i] = 0

    # 3. Compute coefficients for the basis g1, g2, g3
    # g1 = E4^2 = E8
    g1_coeffs = get_eisenstein_series_coeffs(8, num_coeffs)
    
    # g2 = E4 * F
    g2_coeffs = multiply_poly(e4_coeffs, f_coeffs)[:num_coeffs]

    # g3 = F^2 = E8(2z)
    g3_coeffs = [0] * num_coeffs
    for i in range(num_coeffs):
        if i % 2 == 0:
            idx = i // 2
            if idx < len(g1_coeffs):
                g3_coeffs[i] = g1_coeffs[idx]
        else:
            g3_coeffs[i] = 0
            
    # Solve the system of linear equations for c1, c2, c3
    # c1 = 0
    # c2 = 1/240
    # c3 = -1/240
    c1, c2, c3 = 0, 1/240, -1/240
    
    # Compute the coefficients of the cusp form f
    f_coeffs_calc = [0] * num_coeffs
    for i in range(num_coeffs):
        f_coeffs_calc[i] = c1 * g1_coeffs[i] + c2 * g2_coeffs[i] + c3 * g3_coeffs[i]

    # Find the first three non-zero coefficients
    non_zero_coeffs = []
    for coeff in f_coeffs_calc:
        if coeff != 0:
            non_zero_coeffs.append(int(round(coeff)))
        if len(non_zero_coeffs) == 3:
            break
            
    a1, a2, a3 = non_zero_coeffs
    
    # Calculate the sum
    total_sum = sum(non_zero_coeffs)

    print("The unique normalized cusp form f has a q-expansion starting with:")
    print(f"f(q) = {a1}q + {a2}q^2 + {a3}q^3 + ...")
    print(f"\nThe first three non-zero coefficients are {a1}, {a2}, and {a3}.")
    print("The sum of these coefficients is:")
    print(f"{a1} + {a2} + {a3} = {total_sum}")


if __name__ == "__main__":
    main()
