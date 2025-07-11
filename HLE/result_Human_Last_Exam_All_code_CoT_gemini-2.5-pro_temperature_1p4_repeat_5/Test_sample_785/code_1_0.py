import sys

def polymul(c1, c2, max_deg):
    """
    Multiplies two polynomials given as coefficient lists.
    Result is truncated at max_deg.
    """
    deg1 = len(c1) - 1
    deg2 = len(c2) - 1
    new_deg = min(max_deg, deg1 + deg2)
    c_out = [0] * (new_deg + 1)
    
    for i in range(len(c1)):
        if c1[i] == 0:
            continue
        for j in range(len(c2)):
            if c2[j] == 0:
                continue
            if i + j <= max_deg:
                c_out[i + j] += c1[i] * c2[j]
    return c_out

def main():
    """
    Calculates the number of orbits by finding the coefficient of x^1000
    in the corresponding generating function.
    """
    N = 1000
    
    # Explain the problem being solved
    print("The number of orbits is the number of ways to choose non-negative integer")
    print("multiplicities for the 7 irreducible representations of S_5 (with dimensions")
    print("1, 1, 4, 4, 5, 5, 6) such that the total dimension is 1000.")
    print("This is the number of solutions to the equation:")
    print("1*n_1 + 1*n_2 + 4*n_3 + 4*n_4 + 5*n_5 + 5*n_6 + 6*n_7 = 1000")
    print("")

    # Generating function is P(x) = 1/((1-x)^2 * (1-x^4)^2 * (1-x^5)^2 * (1-x^6))
    
    # 1. Coeffs for 1/(1-x)^2 = sum_{k=0 to inf} (k+1)x^k
    coeffs_d1 = [k + 1 for k in range(N + 1)]
    
    # 2. Coeffs for 1/(1-x^4)^2 = sum_{k=0 to inf} (k+1)x^(4k)
    coeffs_d4 = [0] * (N + 1)
    for k in range(N // 4 + 1):
        coeffs_d4[4 * k] = k + 1
        
    # 3. Coeffs for 1/(1-x^5)^2 = sum_{k=0 to inf} (k+1)x^(5k)
    coeffs_d5 = [0] * (N + 1)
    for k in range(N // 5 + 1):
        coeffs_d5[5 * k] = k + 1
        
    # 4. Coeffs for 1/(1-x^6) = sum_{k=0 to inf} x^(6k)
    coeffs_d6 = [0] * (N + 1)
    for k in range(N // 6 + 1):
        coeffs_d6[6 * k] = 1

    # Multiply the polynomial series together
    # Product = (Series_d1)^2 * (Series_d4)^2 * (Series_d5)^2 * Series_d6
    # We can model this by multiplying the coefficient lists.
    # Note: Using the same list twice for squared terms.
    
    term1_coeffs = coeffs_d1  # From 1/(1-x)^2
    term2_coeffs = polymul(term1_coeffs, term1_coeffs, N) # Effectively 1/(1-x)^2 from the two dim 1 irreps
    
    term3_coeffs = polymul(coeffs_d4, coeffs_d4, N) # From the two dim 4 irreps
    
    term4_coeffs = polymul(coeffs_d5, coeffs_d5, N) # From the two dim 5 irreps
    
    term5_coeffs = coeffs_d6 # From the one dim 6 irrep
    
    # Combine terms
    res_coeffs = polymul(term2_coeffs, term3_coeffs, N)
    res_coeffs = polymul(res_coeffs, term4_coeffs, N)
    res_coeffs = polymul(res_coeffs, term5_coeffs, N)
    
    # The answer is the coefficient of x^1000
    result = res_coeffs[N]
    
    print(f"The number of solutions is: {result}")
    
    # The required format for the final answer
    # sys.stdout.write(f"\n<<<{result}>>>")


main()